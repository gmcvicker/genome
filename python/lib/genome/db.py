
import os, sys, re
import tables
import genome.seq

import genome.trackstat
from genome.db_adaptor.chrom_adp import ChromosomeAdaptor

class Track(object):
    def __init__(self, name, path, mode="r"):
        self.name = name
        self.path = path
        self.h5f = tables.openFile(path, mode)

    def __enter__(self):
        sys.stderr.write("Track %s opened\n" % self.name)
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.write("Cleaning up track %s\n" % self.name)
        self.h5f.close()
        return False
    
    def get_array(self, chrom):
        """returns an PyTables ArrayNode of continuous data for a particular
        chromosome"""
        node_name = "/" + str(chrom)
        array_node = self.h5f.getNode(node_name)
        return array_node


    def close(self):
        self.h5f.close()

    
    
    def __get_np_slice(self, array_node, start, end):
        """Helper function, gets a numpy array slice corresponding
        to the provided coordinates from a PyTable ArrayNode"""
        
        if start > end:
            raise ValueError("start (%d) must be <= end (%d)")

        if (start is None) and (end is None):
            return array_node[:]

        if start < 1:
            raise ValueError("start must be >= 1")
        
        if start is None:
            start_idx = 0
        else:
            start_idx = start-1
            
        if end is None:
            end_idx = array_node.shape[0]
        else:
            if end > array_node.shape[0]:
                raise ValueError("end (%d) is greater than chromosome "
                                 "length (%d)" % (end, array_node.shape[0]))
            end_idx = end
        
        return array_node[start_idx:end_idx]


    def get_nparray(self, chrom, start=None, end=None):
        array_node = self.get_array(chrom)
        return self.__get_np_slice(array_node, start, end)


    def get_seq_str(self, chrom, start=None, end=None):
        vals = self.get_nparray(chrom, start, end)
        return genome.seq.from_nparray(vals)




class GenomeDB():
    def __init__(self, path=None):
        if path is None:
            if 'GENOME_DB' in os.environ:
                path = os.environ['GENOME_DB']
            else:
                raise ValueError("No database path specified and "
                                 "GENOME_DB environment variable not set")

        if not path.endswith("/"):
            # add trailing '/' if non exists
            path = path + "/"

        if not os.path.exists(path):
            raise ValueError("database path does not exist: %s" % path)

        if not os.path.isdir(path):
            raise ValueError("database path is not directory: %s" % path)
                
        self.path = path
        
    
    def __enter__(self):
        sys.stderr.write("GenomeDB connection established\n")
        return self


    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.write("Cleaning up GenomeDB\n")
        return False


    def get_track_path(self, track_name):
        if track_name.startswith("/"):
            track_name = track_name[1:]
        if track_name.endswith(".h5"):
            track_name = track_name[:-3]
        track_path = self.path + track_name + ".h5"

        return track_path

        
    def has_track(self, track_name):
        track_path = self.get_track_path(track_name)

        if os.path.exists(track_path):
            return True

        return False

        
    def open_track(self, track_name, mode="r"):
        track_path = self.get_track_path(track_name)

        if not os.path.exists(track_path):
            raise ValueError("track %s does not exist" % track_name)
        
        return Track(track_name, track_path, mode)

    
    def get_track_stat(self, track):
        """Returns statistics for entire track"""        
        return genome.trackstat.get_stats(self, track)

        



    def list_tracks(self, subdir=None):
        """returns a list of existing track names"""
        if subdir is None:
            path = self.path
        else:
            path = self.path + "/" + subdir + "/"
        
        filenames = os.listdir(path)

        track_names = []

        for filename in filenames:
            if filename.startswith("."):
                continue
            if os.path.isdir(path + filename):
                # recursively add tracks to this one
                
                new_subdir = subdir + "/" + filename if subdir else filename
                track_names.extend(self.list_tracks(subdir=new_subdir))
            else:
                if filename.endswith(".h5"):
                    track_name = filename[:-3]
                    if subdir:
                        track_names.append(subdir + "/" + track_name)
                    else:
                        track_names.append(track_name)
        return track_names


    def create_track(self, track_name):
        track_path = self.get_track_path(track_name)

        if os.path.exists(track_path):
            raise ValueError("track %s already exists" % track_path)

        # create parent directories as needed
        dir_names = track_name.split("/")[:-1]
        base_dir = self.path
        for dir_name in dir_names:
            base_dir = base_dir + "/" + dir_name
            if not os.path.exists(base_dir):
                os.mkdir(base_dir)

        return Track(name=track_name, path=track_path, mode="w")


    def get_chromosome_dict(self):
        """Convenience function to retrieve chromosomes in the database,
        returned in a dictionary keyed on chromosome name"""
        chrom_track = self.open_track("chromosome")
        chrom_adp = ChromosomeAdaptor(chrom_track)
        
        return chrom_adp.fetch_name_dict()


    def get_chromosomes(self, get_rand=False, get_auto=True, get_sex=True,
                        get_x=True, get_y=False, get_hap=False,
                        get_mito=False):
        """Convenience function to return chromosomes in the database"""

        chrom_track = self.open_track("chromosome")
        chrom_adp = ChromosomeAdaptor(chrom_track)
        
        return chrom_adp.fetch_by_type(get_rand=get_rand, get_auto=get_auto,
                                       get_sex=get_sex, get_x=get_x,
                                       get_y=get_y, get_hap=get_hap,
                                       get_mito=get_mito)
    


    def get_chromosomes_from_args(self, args):
        """Convenience function, returns a list of chromosome objects
        that are obtained from parsing command line args that may be
        in any one of the following forms 'chr2' or '3', or '1-22'"""

        if len(args) < 1:
            raise ValueError("expected at least one value, got 0")
        
        chrom_name_dict = self.get_chromosome_dict()

        chrom_id_dict = dict([(str(c.idnum), c) for c \
                              in chrom_name_dict.values()])

        chrom_list = []
        for arg in args:
            words = arg.split("-")
            if len(words) == 2:
                # try parsing argument as a range of chromosome
                # ids like 1-22
                try:
                    start = int(words[0])
                    end = int(words[1])

                    if start <= end:
                        vals = [str(x) for x in range(start, end+1)]
                    else:
                        vals = [arg]
                except:
                    vals = [arg]
            else:
                vals = [arg]

            for val in vals:
                if val in chrom_name_dict:
                    chrom_list.append(chrom_name_dict[val])
                elif val in chrom_id_dict:
                    chrom_list.append(chrom_id_dict[val])
                else:
                    raise ValueError("unknown chromosome %s" % val)

        return chrom_list


    def get_chromosome(self, name):
        """Convenience function, retrieves a single chromosome by name"""
        chrom_track = self.open_track("chromosome")
        chrom_adp = ChromosomeAdaptor(chrom_track)
        return chrom_adp.fetch_by_name(name).next()

