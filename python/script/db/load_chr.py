import sys
import tables

import genome.db

# from genome import config
from util import txtfile

class Chromosome(tables.IsDescription):
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


def load_txt_file(filename, chrom_table):
    f = open(filename)
    chrom = chrom_table.row

    for row in txtfile.read_rows(f):
        chrom['idnum'] = int(row["ID"])
        chrom['name'] = row['name']
        chrom['length'] = int(row['length'])
        chrom['is_auto'] = row['is_auto'] == "1"
        chrom['is_sex'] = row['is_sex'] == "1"
        chrom['is_rand'] = row['is_rand'] == "1"
        chrom['is_hap'] = row['is_hap'] == "1"
        chrom['is_mito'] = row['is_mito'] == "1"
        chrom['is_y'] = row['is_y'] == "1"
        chrom['is_x'] = row['is_x'] == "1"
        
        chrom.append()

    chrom_table.flush()    
    f.close()


def main():
    if not len(sys.argv) == 2:
        sys.stderr.write("usage: %s <chrom_txt_file>\n" % sys.argv[0])
        exit(2)

    gdb = genome.db.GenomeDB()
    
    track = gdb.create_track("chromosome")

    chrom_table = track.h5f.createTable("/", 'chromosome', Chromosome,
                                        "chromosomes")

    load_txt_file(sys.argv[1], chrom_table)

    track.h5f.flush()
    track.close()


main()
