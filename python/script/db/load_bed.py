import sys
import tables
import gzip
import argparse

import genome.db
import genome.coord
from genome.coord import CoordError


class Feature(tables.IsDescription):
    start = tables.Int32Col()
    end   = tables.Int32Col()
    strand = tables.Int8Col()
    score = tables.Int16Col()
    name = tables.StringCol(32)



def load_bed_file(filename, chrom_dict, chrom_tab_dict, default_name=".",
                  start_offset=1):
    
    if filename.endswith(".gz"):
        f = gzip.open(filename, "rb")
    else:
        f = open(filename, "r")

    # load the tables with features
    count = 0
    for line in f:
        words = line.rstrip().split()
        
        chrom_name = words[0]
        if chrom_name in chrom_dict:
            chrom = chrom_dict[chrom_name]
        else:
            raise CoordError("unknown chromosome '%s'", chrom_name)

        if chrom_name in chrom_tab_dict:
            chrom_tab = chrom_tab_dict[chrom_name]
        else:
            sys.stderr.write("WARNING: unknown chromosome %s\n" % chrom.name)
            continue
    
        feat = chrom_tab.row

        start = int(words[1]) + start_offset
        end   = int(words[2])

        if len(words) > 5:
            strand = genome.coord.parse_strand(words[5])
        else:
            strand = 0

        if start < 0 or end > chrom.length:
            sys.stderr.write("WARNING: coordinate %d-%d is outside of "
                             "chromosome range 1-%d\n" % (start, end,
                                                        chrom.length))
            continue
        
        if start > end:
            sys.stderr.write("WARNING: start (%d) must be less "
                             "than end (%d)\n" % (start, end))
            continue
        
        if len(words) > 3:
            name = words[3]
            if len(name) > 32:
                sys.stderr.write("WARNING: truncated long name '%s' to '%s'" %
                                 (name, name[0:32]))
                name = name[0:32]
            elif name == ".":
                # use specified name instead
                name = default_name
        else:
            name = default_name

        if len(words) > 5:
            score = int(words[4])
        else:
            score = 0

        # add row to table
        feat['start'] = start
        feat['end'] = end
        feat['strand'] = strand
        feat['name'] = name
        feat['score'] = score    
        feat.append()

        count += 1

    sys.stdout.write("stored %d features\n" % count)
    
    f.close()



def create_track(gdb, chrom_dict, track_name):
    track = gdb.create_track(track_name)
    
    # create separate tables for each chromosome
    chrom_table_dict = {}
    for chrom_name in chrom_dict.keys():
        desc = track.name + " features for " + chrom_name
        chrom_tab = track.h5f.createTable("/", chrom_name, Feature, desc)
        chrom_table_dict[chrom_name] = chrom_tab

    return track, chrom_table_dict
    

def main():
    parser = argparse.ArgumentParser()    
    parser.add_argument("--start_offset", type=int, default=1,
                        help="value to add to start coordinate")
    parser.add_argument("track_name", help="name of track to store coords in")

    parser.add_argument("filenames", help="BED-like files to read coords from",
                        nargs='+')

    options = parser.parse_args()

    gdb = genome.db.GenomeDB()

    chrom_dict = gdb.get_chromosome_dict()
    track, chrom_tab_dict = create_track(gdb, chrom_dict, options.track_name)
    
    for filename in options.filenames:
        load_bed_file(filename, chrom_dict, chrom_tab_dict,
                      start_offset=options.start_offset)

    # flush tables
    for tab in chrom_tab_dict.values():
        tab.flush()

    track.h5f.flush()
    track.close()



if __name__ == "__main__":
    main()
