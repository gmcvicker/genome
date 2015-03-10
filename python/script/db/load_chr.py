import sys
import gzip
import re
import tables
import argparse

import genome.db
import genome.chrom

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


def load_chromosomes(chrom_list, chrom_table):
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



def chrom_key(chrom):
    """Returns a key for sorting chromosomes based on their name"""
    m = re.match(r"^chr(\d+)", chrom.name)    
    if m:
        # make sure autosomes are sorted numerically by padding with
        # leading 0s
        num = m.groups()[0]

        if len(num) < 3:
            name = ("0" * (3-len(num))) + num
        else:
            name = num
    else:
        # otherwise just sort lexigraphically
        name = chrom.name

    # first take non-haplo, non-rand, non-sex chromosomes, then
    # sort by name
    return (chrom.is_hap, chrom.is_rand, chrom.is_mito, chrom.is_sex, name)



def parse_chromosomes(filename):
    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)

    chrom_list = []
    
    for line in f:
        words = line.rstrip().split()

        if len(words) < 2:
            raise ValueError("expected at least two columns per line\n")
        
        chrom = genome.chrom.Chromosome(name=words[0], length=words[1])
        chrom_list.append(chrom)

        lc_name = chrom.name.lower()

        # determine whether this is autosome, sex or mitochondrial chrom
        if re.match('^chr(\d+)', lc_name):
            chrom.is_auto=True
        elif re.match("^chr[W-Zw-z]", lc_name):
            chrom.is_sex = True
        elif lc_name.startswith("chrm"):
            chrom.is_mito = True
        elif lc_name.startswith("chrun") or lc_name.startswith("chrur"):
            chrom.is_rand = True
        else:
            sys.stderr.write("WARNING: could not determine chromosome type "
                             "(autosome, sex, mitochondrial) from name "
                             "'%s'. Assuming 'random'\n" % chrom.name)
            chrom.is_rand = True

        if "rand" in chrom.name:
            # random chromosome
            chrom.is_rand = True

        if "hap" in chrom.name:
            # alt haplotype chromosome
            chrom.is_hap = True

    chrom_list.sort(key=chrom_key)

    idnum = 1
    for chrom in chrom_list:
        chrom.idnum = idnum
        idnum += 1

        sys.stderr.write("%s\n" % chrom.name)

    f.close()

    return chrom_list



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", default=None,
                        help="assembly to create chromosome table for"
                        " (e.g. hg19)")

    parser.add_argument("chrom_file",
                        help="path to UCSC chromInfo.txt.gz file")

    return parser.parse_args()
    


def main():
    args = parse_args()

    sys.stderr.write("parsing input file %s\n" % args.chrom_file)
    chrom_list = parse_chromosomes(args.chrom_file)

    gdb = genome.db.GenomeDB(assembly=args.assembly)

    sys.stderr.write("creating chromosome table for %s\n" %
                     gdb.assembly)
    track = gdb.create_track("chromosome")

    chrom_table = track.h5f.createTable("/", 'chromosome', ChromDesc,
                                        "chromosomes")

    sys.stderr.write("storing chromosomes in table\n")
    load_chromosomes(chrom_list, chrom_table)

    track.h5f.flush()
    track.close()


main()
