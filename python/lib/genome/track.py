import sys
import os
import tables
import numpy as np
import datetime

import genome.seq
import genome.trackstat


class ChromDesc(tables.IsDescription):
    idnum = tables.Int32Col()
    name = tables.StringCol(32)
    length = tables.Int32Col()
    is_sex = tables.BoolCol(dflt=False)
    is_auto = tables.BoolCol(dflt=True)
    is_rand = tables.BoolCol(dflt=False)
    is_hap = tables.BoolCol(dflt=False)
    is_mito = tables.BoolCol(dflt=False)
    is_y = tables.BoolCol(dflt=False)
    is_x = tables.BoolCol(dflt=False)


class MetaDataDesc(tables.IsDescription):
    key = tables.StringCol(32)
    val = tables.StringCol(1024)



class Track(object):
    """This class represents a data track and is an abstraction over a
    single HDF5 file. It allows for easy retrieval of a numpy array or
    DNA sequence string for a set of chromosomal coordinates. It can
    also be used to obtain meta information about the track.

    In theory this class could be exteded to work as an abstraction over
    several different kinds of file types (e.g. bigWig or bam),
    although I'm not certain this is a good idea.
    """
    
    def __init__(self, path, mode="r", chromosomes=None,
                 assembly=None, name=None, desc=None):
        self.path = path

        if mode == "r":
            self.h5f = tables.open_file(path, "r")
        elif mode == "w":
            # creating a new file
            if os.path.exists(path):
                raise IOError("file %s already exists, "
                              "please remove before opening in 'w' mode" %
                              path)

            if chromosomes is None:
                raise ValueError("must provide list of chromosomes "
                                 "when creating new h5f track")

            if assembly is None:
                raise ValueError("must provide name of assembly "
                                 "when creating new h5f track")

            if name is None:
                raise ValueError("must provide name of track "
                                 "when creating new h5f track")

            if desc is None:
                raise ValueError("must provide description of track "
                                 "when creating new h5f track")
            
            self.h5f = tables.open_file(path, mode)            
            self.load_chromosomes(chromosomes)
            
            #
            # load meta information into meta table
            #            
            meta_dict = {'assembly' : assembly,
                         'name' : name,
                         'desc' : desc}
            
            if "USER" in os.environ:
                meta_dict['user'] = os.environ['USER']
                            
            meta_dict['command'] = " ".join(sys.argv)
            meta_dict['date'] = datetime.datetime.now().strftime("%b %d, %Y, %H:%M")

            self.load_meta(meta_dict)
        
        # read meta information from meta table
        meta_dict = self.get_meta_dict()

        self.meta = meta_dict

        if 'assembly' in meta_dict:
            self.assembly = meta_dict['assembly']
        else:
            sys.stderr.write("WARNING: assembly not specified in track\n")
            self.assembly = None

        if 'name' in meta_dict:
            self.name = meta_dict['name']
        else:
            sys.stderr.write("WARNING: name not specified in track\n")
            self.name = None

        if 'desc' in meta_dict:
            self.desc = meta_dict['desc']
        else:
            sys.stderr.write("WARNING: desc not specified in track\n")
            self.desc = None
                
        self._missing_chrom = set([])


        
    def __enter__(self):
        sys.stderr.write("Track %s opened\n" % self.name)
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.write("Cleaning up track %s\n" % self.name)
        self.h5f.close()
        return False

    
    def load_chromosomes(self, chrom_list):
        """Creates a table containing list of chromosomes and 
        their lengths"""        
        chrom_table = self.h5f.create_table("/", 'chromosome', ChromDesc,
                                            "chromosomes")

        row = chrom_table.row

        for chrom in chrom_list:
            row['idnum'] = chrom.idnum
            row['name'] = chrom.name
            row['length'] = chrom.length
            row['is_auto'] = chrom.is_auto
            row['is_sex'] = chrom.is_sex
            row['is_rand'] = chrom.is_rand
            row['is_hap'] = chrom.is_hap
            row['is_mito'] = chrom.is_mito
            row['is_y'] = chrom.is_y
            row['is_x'] = chrom.is_x                
            row.append()

        chrom_table.flush()



    def load_meta(self, meta_dict):
        """creates a table containing meta data key/value pairs"""
        meta_table = self.h5f.create_table("/", 'meta', MetaDataDesc,
                                          "meta data")

        row = meta_table.row

        for k, v in meta_dict.items():
            row['key'] = k
            row['val'] = str(v)
            row.append()

        meta_table.flush()


    def get_meta_dict(self):
        """reads meta key/value pairs from table"""
        meta_dict = {}

        for row in self.h5f.root.meta:
            key = row['key'].decode("utf-8")
            val = row['val'].decode("utf-8")
            meta_dict[key] = val
            sys.stderr.write("%s => %s\n" % (key, val))


        return meta_dict


    def add_meta_val(self, key, val):
        """Adds a value to the meta table. currently raises error 
        if table not opened in append mode, or if key is already in
        table"""
        
        meta_dict = get_meta_dict()

        if key in meta_dict:
            raise ValueError("key %s is already in meta table" % key)

        table = self.h5f.root.meta
        row = table.row
        row['key'] = key
        row['val'] = str(val)
        row.append()

        table.flush()

    
    def get_chromosome_dict(self):
        """Returns a dictionary of all chromosomes in the database,
        keyed on chromosome name"""
        chrom_list = self.get_all_chromosomes()
        chrom_dict = {}
        for chrom in chrom_list:
            chrom_dict[chrom.name] = chrom

        return chrom_dict
        


    def get_chromosomes(self, get_rand=False, get_auto=True, get_sex=True,
                        get_x=True, get_y=False, get_hap=False,
                        get_mito=False):
        """Returns a filtered list of Chromosomes from the
        database. Optional flags specify the subset of Chromosomes
        that are returned. By default the 22 autosomes and chrX are
        retrieved (but chrY, the mitochondrial chromosome, alternate
        haplotypes, and 'random' chromosomes are not)"""
        chrom_list = []
        chrom = None

        if not "chromosome" in self.h5f.root:
            raise ValueError("need to load chromosome table by calling "
                             "load_chromosomes first")


        flags = {'is_rand' : get_rand,
                 'is_auto' : get_auto,
                 'is_sex' : get_sex,
                 'is_x' : get_x,
                 'is_y' : get_y,
                 'is_hap' : get_hap,
                 'is_mito' :  get_mito}

        conditions = []
        for k, v in list(flags.items()):
            if not v:
                # don't get this chromosome type
                conditions.append("(%s == False)" % k)

        row_iter = None
        if len(conditions) > 0:
            query = " & ".join(conditions)
            row_iter = self.h5f.root.chromosome.where(query)
        else:
            row_iter = self.h5f.root.chromosome

        for row in row_iter:
            chrom = Chromosome(idnum=row['idnum'],
                               name=row['name'],
                               length=row['length'],
                               is_auto=row['is_auto'],
                               is_hap=row['is_hap'],
                               is_mito=row['is_mito'],
                               is_sex=row['is_sex'],
                               is_x=row['is_x'],
                               is_y=row['is_y'])

            chrom_list.append(chrom)

        return chrom_list
    


    
    def get_all_chromosomes(self):
        """Returns an unfiltered list of all of the chromosomes in the
        database"""
        chrom_list = []
        chrom = None

        if not "chromosome" in self.h5f.root:
            raise ValueError("need to load chromosome table by calling "
                             "load_chromosomes first")

        for row in self.h5f.root.chromosome:
            chrom = Chromosome(idnum=row['idnum'],
                               name=row['name'],
                               length=row['length'],
                               is_auto=row['is_auto'],
                               is_hap=row['is_hap'],
                               is_mito=row['is_mito'],
                               is_sex=row['is_sex'],
                               is_x=row['is_x'],
                               is_y=row['is_y'])

            chrom_list.append(chrom)

        return chrom_list
        

        


    def get_chromosomes_from_args(self, args):
        """Convenience function, returns a list of chromosome objects
        that are obtained from parsing command line args that may be
        in any one of the following forms 'chr2' or '3', or '1-22'"""
        if isinstance(args, str):
            # if we were given a string, put it in a list
            args = [args]
        
        if len(args) < 1:
            raise ValueError("expected at least one value, got 0")
        
        chrom_name_dict = self.get_chromosome_dict()

        chrom_id_dict = dict([(str(c.idnum), c) for c \
                              in list(chrom_name_dict.values())])

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
        """Retrieves a single chromosome by name"""
        chrom_dict = self.get_chromosome_dict()
        return chrom_dict[name]
        
    
    def has_chromosome(self, chrom):
        """Returns True if this track contains a particular chromosome,
        False otherwise"""
        chrom_dict = self.get_chromosome_dict()
        return chrom in self.chromosome_dict()

    
    
    def get_array(self, chrom):
        """returns an PyTables ArrayNode for a particular chromosome"""
        node_name = "/" + str(chrom)

        if node_name in self.h5f:
            array_node = self.h5f.get_node(node_name)
        else:
            if str(chrom) not in self._missing_chrom:
                sys.stderr.write("WARNING: track '%s' is missing "
                                 "chromosome '%s'\n" % 
                                 (self.name, str(chrom)))
                self._missing_chrom.add(str(chrom))
            return None
                
        return array_node

    
    def get_val(self, chrom, pos):
        """returns the value of the track at the specified "
        genomic position"""

        array = self.get_array(chrom)

        if array:
            return array[pos-1]

        return np.nan

    
    def close(self):
        """Closes this track by closing the underlying HDF5 file"""
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
        """Returns a numpy array of data for the specified chromosome
        or chromosomal region"""
        array = self.get_array(chrom)

        if array is None:
            # default to array of nan
            if hasattr(chrom, "length"):
                array = np.empty(chrom.length, dtype=np.float32)
            else:
                raise ValueError("cannot create array for missing chromosome "
                                 "of unknown length for track '%s'" % self.name)

            array[:] = np.nan

        return self.__get_np_slice(array, start, end)

        

    def get_seq_str(self, chrom, start=None, end=None):
        """Returns a string of sequence of the specified chromosome
        or chromosomal region. It only makes sense to call this function
        for tracks represent sequence data as 8-bit integers that can be
        converted to printable characters."""
        vals = self.get_nparray(chrom, start, end)
        return genome.seq.from_nparray(vals)


    
    def get_stats(self, track):
        """Returns a TrackStat object containing statistics for
        this track (mean, max, sum, etc.)"""
        return genome.trackstat.get_stats(track)


    def set_stats(self, track):
        """computes and sets track statistics for this track object"""
        genome.trackstat.set_stats(self, track)




