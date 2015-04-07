import sys
import gzip

import argparse

from genome.transcript import Transcript
from genome.coord import Coord
import genome.db

GENE_TYPE_KNOWN = 1
GENE_TYPE_ENS = 2
GENE_TYPE_REF = 3
GENE_TYPE_CCDS = 4
GENE_TYPE_FLYBASE = 5


def merge_adjacent_exons(exons):
    new_exons = []

    cur_exon = exons[0]
    new_exons.append(cur_exon)

    for exon in exons[1:]:
        if exon.start == cur_exon.end+1:
            cur_exon.end = exon.end
        else:
            cur_exon = exon
            new_exons.append(cur_exon)

    return new_exons


def write_header():
    print("\t".join(["ID", "NAME", "CHROM", "START", "END", "STRAND",
                     "EXON.STARTS", "EXON.ENDS",
                     "CDS.START", "CDS.END"]))


def parse_args():
    parser = argparse.ArgumentParser(description="parses several types "
                                     "of UCSC gene formats, outputs a "
                                     "common format for all of them")

    parser.add_argument("--assembly", metavar="ASSEMBLY", default=None,
                        help="name of assembly (e.g. hg18)")

    parser.add_argument("gene_file", metavar="GENE_FILE",
                        help="path to ucsc gene file")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    path = args.gene_file
    filename = path.split("/")[-1]

    if filename.startswith("known"):
        gene_type = GENE_TYPE_KNOWN
    elif filename.startswith("ens"):
        gene_type = GENE_TYPE_ENS
    elif filename.startswith("ref"):
        gene_type = GENE_TYPE_REF
    elif filename.startswith("ccds"):
        gene_type = GENE_TYPE_CCDS
    elif filename.startswith("flyBase"):
        gene_type = GENE_TYPE_FLYBASE
    else:
        raise ValueError("cannot discern input gene type from filename")

    if path.endswith(".gz"):
        f = gzip.open(path)
    else:
        f = open(path)

    gdb = genome.db.GenomeDB(assembly=args.assembly)

    chrom_dict = gdb.get_chromosome_dict()

    tr_id = 0

    write_header()

    skipped_chrom = []

    for l in f:
        tr_id += 1
        cols = l.strip().split("\t")

        if gene_type in (GENE_TYPE_ENS, GENE_TYPE_REF, GENE_TYPE_CCDS,
                         GENE_TYPE_FLYBASE):
            # skip first column, which is "bin"
            cols = cols[1:]

        name = cols[0]

        chrom_name = cols[1]
        if not chrom_name in chrom_dict:
            if not chrom_name in skipped_chrom:
                sys.stderr.write("skipping unknown chromosome: %s\n"
                                 % chrom_name)
                skipped_chrom.append(chrom_name)
            continue
        chrom = chrom_dict[chrom_name]

        if cols[2] == '+':
            strand = 1
        elif cols[2] == '-':
            strand = -1
        else:
            sys.stderr.write("WARNING: unknown strand '%s' for gene "
                             "%s--skipping\n" % (cols[2], name))
            # raise ValueError("invalid strand: %s" % cols[2])

        cds_start = int(cols[5])+1
        cds_end = int(cols[6])

        if(cds_start == cds_end+1):
            # indicates non-coding transcript
            cds_start = None
            cds_end = None

        start_strs = cols[8].split(",")
        end_strs = cols[9].split(",")
        start_strs.pop()
        end_strs.pop()
        exon_starts = [int(x)+1 for x in start_strs]
        exon_ends = [int(x) for x in end_strs]

        exons = []
        for start, end in zip(exon_starts, exon_ends):
            exon = Coord(chrom, start, end, strand)
            exons.append(exon)

        # merge UCSC exons that are immediately adjacent
        # i.e. 0 bp introns
        exons = merge_adjacent_exons(exons)

        if strand == -1:
            exons.reverse()

        tr = Transcript(name=name, exons=exons,
                        cds_start=cds_start, cds_end=cds_end, idnum=tr_id)

        print(str(tr))

    f.close()


main()
