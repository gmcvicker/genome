import tables
import pysam
import sys
import genome.db


MIN_READ_LEN = 1
MAX_READ_LEN = 500

MIN_MAP_QUAL = 10

class Feature(tables.IsDescription):
    start = tables.Int32Col()
    end = tables.Int32Col()
    strand = tables.Int8Col()
    score = tables.Int8Col()


def load_bam_reads(chrom, chrom_tab, bam_filename):
    samfile = pysam.Samfile(bam_filename, "rb")

    # create separate table for each chromosome
    count = 0

    for read in samfile.fetch(reference=chrom.name,
                              start=1, end=chrom.length):
        count += 1
        if count > 100000:
            sys.stderr.write(".")
            count = 0
        
        # reads appear twice, once for each side, only want to consider once
        if not (read.is_read1 or read.is_read2):
            continue
        if read.is_read2:
            continue

        # require that both sides are uniquely mapped
        if read.is_unmapped or read.mate_is_unmapped:
            continue

        if read.is_reverse == read.mate_is_reverse:
            # reads mapped to same strand...
            continue

        if read.mapq < MIN_MAP_QUAL:
            # read has poor mapping quality
            continue

        # remember pysam pos starts at 0, not 1
        if read.is_reverse:
            isize = -read.isize
            start = read.mpos + 1
        else:
            isize = read.isize
            start = read.pos + 1
        
        if isize < MIN_READ_LEN or isize > MAX_READ_LEN:
            continue

        if read.is_reverse:
            strand = -1
        else:
            strand = 1

        end = start + isize - 1
        
        if start < 0 or end > chrom.length:
            sys.stderr.write("skipping read %d-%d, outside of chromosome "
                             "range 1-%d" % (start, end, chrom.length))
            continue
        if start > end:
            raise CoordError("start (%d) must be less than end (%d)" %
                             (start, end))


        feat = chrom_tab.row
        feat['start'] = start
        feat['end'] = end
        feat['strand'] = strand
        feat['score'] = read.mapq

        feat.append()



def main():
    if len(sys.argv) < 3:
        sys.stderr.write("usage: %s <track_name> <bamfile1> "
                         "[<bamfile2> ...])\n" % sys.argv[0])
        exit(2)

    track_name = sys.argv[1]
    filenames = sys.argv[2:]

    gdb = genome.db.GenomeDB()

    track = gdb.create_track(track_name)

    for chrom in gdb.get_chromosomes():        
        # create a new feature table for this chromosome        
        sys.stderr.write("%s\n" % chrom.name)
        desc = track.name + " reads for " + chrom.name
        chrom_tab = track.h5f.createTable("/", chrom.name, Feature, desc)

        # load aligned read coordinates for each BAM file
        for filename in filenames:
            sys.stderr.write("  %s\n  " % filename)
            load_bam_reads(chrom, chrom_tab, filename)
            sys.stderr.write("\n")
            
        
        chrom_tab.flush()

    track.close()
            
main()
