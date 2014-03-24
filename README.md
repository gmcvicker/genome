# Introduction

This repo contains a Python library, *genome*, for retrieval and storage of genomics data in 
[HDF5](http://www.hdfgroup.org/HDF5/) format. It is
a lightweight wrapper over [PyTables](http://www.pytables.org/moin).

This repo contains the *genome* library and a set of utility scripts for reading and 
writing [HDF5](http://www.hdfgroup.org/HDF5/) files. 

The purpose of this software is to:
* Provide a simple abstraction over PyTables for retrieval of genomics data
* Make it simple to convert other genomics data formats to HDF5

This library was inspired by and is loosely based on the 
[Genomedata](http://noble.gs.washington.edu/proj/genomedata/) software
by [Michael Hoffman](http://noble.gs.washington.edu/~mmh1/).

This document summarizes how to setup this software system, how to retrieve data from HDF5 files, 
and how to convert other file formats into HDF5. If you have any questions or run into any 
difficulty, please don't hesitate to contact me!

### HDF5 and PyTables

HDF5 is a binary format designed for rapid access to huge amounts of data. 
It is widely-adopted and platform independent. Software libraries to access HDF5 file are available in 
C, C++, Java, Python, and Fortran90 (and possibly others). PyTables is a very nice high-level Python 
interface for reading and writing HDF5 files. I have written a python library and a number of scripts 
that use PyTables to create and read HDF5 files.

Advantages to using HDF5
* Efficient storage of very large datasets
* Very fast data retrieval
* Fast random access to compressed files
* Nice, simple Python interface (PyTables) that integrates with NumPy
* APIs in several languages: C, C++, Java, Fortran
* Binary files work on big- and little-endian machines

Disadvantages:
* Files are not human readable
* PyTables does not support write-concurrency (only one process at a time can safely write to an HDF5 file)
* Unlike MySQL and PostgreSQL, The free version of PyTables does not support indexing  (the commercial version does)
* Using PyTables introduces several software dependencies: PyTables, NumPy, HDF5, zlib


# Setup

This repo has three components: a small python library *genome*, a small C 
library *libgenome* and set of python scripts for data import and manipulation.

### Dependencies

*genome* depends on PyTables, which depends on HDF5(http://www.hdfgroup.org/HDF5/), [NumPy](http://www.numpy.org/) and several other packages. 
See the detailed [PyTables installation guide](http://pytables.github.io/usersguide/installation.html). [This 
guide](http://assorted-experience.blogspot.com/2011/12/mac-os-x-install-pytables-and-h5py.html) for installing
PyTables under Mac OSX.

*genome* additionally depends on:
* [Cython](http://cython.org/)
* [ArgParse](https://code.google.com/p/argparse/) (for Python2.6 or older)
* [pysam](https://code.google.com/p/pysam/) (for reading BAM files only)
* [Scons](http://www.scons.org/) (for compiling, would like to replace with standard build scripts)


### Getting the source code

The first step is to obtain the source code from git. Let's assume you want to clone this repository 
into a directory named ~/src:

		mkdir ~/src
		cd ~/src
		git clone https://github.com/gmcvicker/genome

If at a later point you want to retrieve the latest code updates you can do the following:

		cd ~/src/genome
		git fetch
		git merge origin/master


### Setting environment variables

To use the code you need to update set several shell environment variables. 
I recommend setting these variables in your ~/.bashrc or ~/.profile files using 
something like the following:

    # update your python path by adding $HOME/src/genome/python/lib to the end
    # this tells python where to find the genome library 
    export PYTHONPATH=$PYTHONPATH:$HOME/src/genome/python/lib

    # update LD_LIBRARY_PATH by adding $HOME/src/genome/c/lib/ 
    # this is so the libgenome C library can be dynamically loaded by scripts that need it
    export LD_LIBRARY_PATH=$HOME/src/genome/c/lib/:$LD_LIBRARY_PATH
    
    # set some flags that are needed for compiling the C library
    export CFLAGS=-I$HOME/src/genome/c/lib/
    export LDFLAGS=-L$HOME/src/genome/c/lib/
    
    # specify the location of 'database' where the HDF5 files that you want to use are
    export GENOME_DB=/data/share/genome_db/
    
    # scripts will use this genome assembly by default 
    export GENOME_ASSEMBLY=hg19

Make sure that you remember to set these variables after adding them to your .bashrc for the first time. 
You can login again, or do source ~/.bashrc

Remember that when you submit jobs to a compute cluster (e.g. using SGE's qsub), they run in "batch mode" 
and may not execute your ~/.bashrc. To ensure that your jobs have the correct environment variables set, 
you should be able to pass a flag to your cluster submission command (e.g. the -V flag to qsub).


### Compiling the C library and building python extensions

I use a build system called [SCons](http://www.scons.org/) to compile my code. To compile the libgenome 
C library type the following:

    cd ~/src/genome/c
    scons


The genome library uses [Cython](http://cython.org/) to create C extensions for python. 
To build these extensions do the following:

    # build C extensions that are part of genome library
    cd ~/src/genome/python/lib/genome
    python setup.py build_ext --inplace

    # build C extensions that are needed by import scripts
    cd ~/src/genome/python/script/db/
    python setup.py build_ext --inplace

This may give some warnings, which are OK. If you get errors that indicate the build failed, 
this might be because the LD_LIBRARY_PATH, CFLAGS, LDFLAGS environment variables are set incorrectly, 
or because the genome C library is not properly compiled.



# Importing Data

A number of scripts can be used to load data into the HDF5 files. Most of these scripts will 
take a name for a new track. The data will then be written to file with the name of the track 
(and a .h5 postfix) under the directory pointed to by the GENOME_DB environment variable.

Subdirectories are automatically created for tracks. For example, if you create a new track called 
dnase/new_dnase_track this will create a new HDF5 file named 
$GENOME_DB/$GENOME_ASSEMBLY/dnase/new_dnase_track.h5. These scripts will report an error if you 
try to write over an existing track with the same name. If such a file exists you will need to remove it manually.

The scripts to convert files into HDF5 are located in genome/python/script/db. The 
create_track.py script is the most flexible, and it can take several different 
file formats as input. The following summarizes the scripts can be used to import different types of data.

*Note:* Many of these scripts require the [argparse](http://docs.python.org/2/howto/argparse.html) module. 
argparse is part of python2.7, but if you are using python2.6 or earlier then you will have to 
install argparse on your system.

*Note:* The scripts that read BAM files depend on [pysam](https://code.google.com/p/pysam/). 

### Creating a new database

A new database needs only one track, the chromosome track (which is actually more a table than a track). 
This can be created using the load_chr.py script described below. It may also be desirable to load a genome
sequence track. An example of loading a genome sequence from fasta files is given in the
create_track.py section below.


### load_chr.py

Creates a table of chromosomes in the database. This table contains chromosome names, lengths 
and some other information. This is the only table that is required for the libary to work, and it 
only needs to be created once. If you want to create a new database (e.g. for 
another species or on a different machine) you will have to use this script to create a new table, 
or copy an existing chromosome.h5 file. The load_chr.py script takes a chromInfo.txt.gz file as
input. These can be downloaded from the UCSC genome browser.


### create_track.py

Takes fasta, wiggle, bedgraph, xb or tab-delimited text files as input. Data are stored in a track containing a 
1D array for each chromosome. The imported data can currently be stored as any of the 
following data types: int8, uint8, int16, float32. Currently create_track.py expects 
fasta, wiggle, bedgraph and text files to be split into separate chromosomes with 
filenames that contain the chromosome name (e.g. chr1.wig.gz, chr2.wig.gz, etc.). 
The provided C program genome/c/program/split_wig_chrs.c can be used to split wiggle files up in this way.

Here is an example of how to load sequence data into the database (in this case for Drosophila melanogaster):

    python create_track.py --assembly dm3 --format fasta --dtype uint8 seq  ~/data/Dmel/ucsc/dm3/seq/chr*.fa.gz

### load_bed.py

Reads features from a BED file and stores them in a HDF5 file. Data imported this way 
are stored in a table with several columns (not as a 1D array).

### load_bam_read_depth.py

Reads BAM or SAM files and stores read depths in a track as a 1D array of unsigned 
16 bit integers (uint16s) for each chromosome. Discards strand information but could easily 
be modified to store forward and reverse strands separately. Requires the pysam python library. 
BAM file must first be sorted and indexed using [samtools](http://samtools.sourceforge.net/).
 
### load_bam_5prime_ends.py

Reads BAM or SAM files and stores counts of 5' ends of reads as a 1D array of unsigned 16 bit 
integers (uint16s) for each chromosome. Forward and reverse strands are stored as separate tracks. 
Requires the pysam python library. BAM file must first be sorted and indexed using 
[samtools](http://samtools.sourceforge.net/). 
Example of use (note any number of bam files can be specified):

    python load_bam_5prime_ends.py --assembly hg19 \
      	/my_tracks/fwd_read_counts \
        /my_tracks/rev_read_counts \
        wgEncodeDataRep1.bam   wgEncodeDataRep2.bam

### load_bam_left_ends.py

Reads BAM or SAM files and stores counts of the left (lower-numbered-coordinate) ends of reads 
as a 1D array of unsigned 16 bit integers (uint16s) for each chromosome. Forward and reverse 
strands are stored in the same tracks. Requires the pysam python library. BAM file must first 
be sorted and indexed using samtools. Example of use:

    python load_bam_left_ends.py --assembly hg19 /my_tracks/read_counts \
        wgEncodeDataRep1.bam wgEncodeDataRep2.bam wgEncodeDataRep3.bam

### load_mnase_mids.py

Estimates dyad positions from single-end or paired-end MNase-seq reads and stores dyad counts as 1D
arrays of unsigned 8 bit integers. Data are read from a BAM file, which must first be sorted and 
indexed using samtools. Requires the pysam python library. 



# Retrieving Data

You should now be ready to obtain data using the genome library. The environment variable 
GENOME_DB tells genome.db where your database of HDF5 files is located. 
Each .h5 file under this directory (or within subdirectories) is a "track" from which data can be 
retrieved. The HDF5 files that I have created reside in /data/share/genome_db/ so you
should make sure that your GENOME_DB environment variable points there.

The following script is an example showig how to retrieve DNase and MNase data from the database. 
Note that this script assumes that your database already contains tracks named 'dnase/dnase_all_combined'
and 'mnase/q10/mnase_mids_combined'. Creation of tracks and importing of data is described below.

    import numpy as np
    import genome.db
    
    # 'connect' to the genome database
    gdb = genome.db.GenomeDB()
    
    # open data 'tracks' for DNase and MNase
    dnase_track = gdb.open_track('dnase/dnase_all_combined')
    mnase_track = gdb.open_track('mnase/q10/mnase_mids_combined')

    # pull out numpy arrays for a region of chromosome 12
    start = 1000000
    end   = 1000100
    dnase_vals = dnase_track.get_nparray('chr12', start, end)
    mnase_vals = mnase_track.get_nparray('chr12', start, end)
    print dnase_vals
    print mnase_vals
    
    # loop over chromosomes 1-22 and X counting the totals on each
    for chrom in gdb.get_chromosomes():
	    # report the length of the chromosome
	    print "%s length is %dbp" % (chrom.name, chrom.length)

      # Count dnase cuts on this chromosome.  Ask for the entire
      # chromosome's values at once.  This uses a lot more memory than
      # getting small regions, but is faster if we want to look at every
      # site
      vals = dnase_track.get_nparray(chrom.name)
      total = np.sum(vals)
    print "DNase cuts: %d" % total

    # now count number of MNase midpoints
	  vals = mnase_track.get_nparray(chrom.name)
	  total = np.sum(vals)
	  print "MNase midpoints: %d" % total
    
    # close the HDF5 files
    mnase_track.close()
    dnase_track.close()

The output from this program looks like this. This program normally only takes a couple of minutes to run.
If it runs very slowly this is probably because somebody else has submitted a lot
of cluster jobs that are taxing the lustre filesystem.

    [3 5 3 2 3 1 2 0 0 1 0 0 0 0 2 2 1 0 2 2 1 1 0 1 0 0 0 0 4 2 4 1 0 0 1 1 1 
     2 0 0 1 1 1 3 1 0 0 0 0 0 1 3 0 0 0 2 1 1 1 0 3 3 2 1 0 1 0 2 3 5 2 1 5 0 
     2 0 0 0 1 0 1 0 1 4 0 3 1 2 5 1 0 1 4 0 2 1 0 0 2 0 1]
    [7  6  2  2  2  4  2  0  0  3  5  2  1  0  2  2  4  0  1  3  4  3 28  1  6
     0  5  2 19  4  1  0  0  0  1  0  0  2  1  1  5  6  1  3  4  1  2  9  7  2
     2  2  3  1  0  0  0  2  0  0  0  0  4  2  2  4  4  1  3 48  0  2  0  1  0
     1  0  0  1  1  1  3  0  0  1  1  1  2  0  1  0  0  0  1  9  1  0  1  0 25
     2]
    chr1 length is 247249719bp
    DNase cuts: 265626609
    MNase midpoints: 301957846
    chr2 length is 242951149bp
    DNase cuts: 236830314
    MNase midpoints: 311127578
    ...

Here is another example. This program retrieves DNA sequence for a genomic region that is 
specified on the command line.

    if len(sys.argv) != 4 and len(sys.argv) != 5:
	    sys.stderr.write("usage: %s <chr> <start> <end> [+|-]\n" % sys.argv[0])
    exit(2)

    # parse command line arguments
    chrom_name = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])

    if len(sys.argv) == 5:
	    strand = genome.coord.parse_strand(sys.argv[4])
    else:
	    strand = 1

    # 'connect' to database
    gdb = genome.db.GenomeDB()

    # get a track which contains DNA sequence
    track = gdb.open_track("seq")

    # retrieve a sequence string from this track
    seq_str = track.get_seq_str(chrom_name, start, end)

    if strand == -1:
	    # reverse complement the sequence if the reverse strand was requested
	    seq_str = genome.seq.revcomp(seq_str)
	    strand_str = '(-)'
    else:
	    strand_str = ''

    # write FASTA-formatted output
    seq_id = "%s:%d-%d%s" % (chrom_name, start, end, strand_str)
    genome.fasta.write_fasta(sys.stdout, seq_id, seq_str)

    # close the HDF5 file
    track.close()

The output from running this program looks like this:

    python get_seq.py chr12 1000000 1010000
    >chr12:1000000-1010000
    CACTCATTTCAATTCTTGACTGGGATTGGAGGATCTTTGAAGCTATCTCACGTGGCTGTT
    GGCAGGAGACCTCAGTTTCTGAATATGTGACTCTCTTCACAGGGTAGTTTGAAGTTTGAG
    TATCCTCACGACATGGCAGCTGGATTCCTGAGAGAGAGTGATTTGAGAGAGTGACGAAGA
    CAGAAGCCCCAGTGTCTTTTAAGAGTTAGTCATGAAAGTACATCATTGGAAGCAAGTCAT
    TAAGTCCAGTCCACACTCAAGGAGAGGAAATTTAGTCTCCCTCTTGAAGGGAGGAGCATT
    AAATAATTTATGAACATACTTTAAAACTACCACAAATAATTGTCAGGAGAATGAATCTTC
    CATGAGCTTCTTATTTCTTTTGTTTGAGATGGAGTCCTGCTCTTTCACCCAGGCAGGAGT
    GCAGTGGCGCAGTCTGGGCTCACTGCAGCCTCCACCTCCCGGGTTCAAGCGATTCTTCTG
    CCTCAGCCTCCTGAGTAGCTGGGATTACAGGCACATGGCCACCATGCCTGGCTAATTTTT
    GTATTTTTAGTAGAGACGGGGTTTCACCATGTTGGCCAAACTGGTCTCAAACTCCTGACC
    ...



# Other scripts

Here are some other scripts that may be quite useful. They are also located in genome/python/script/db. 

#### list_tracks.py
Print a list of tracks that are in the database

#### list_chromosomes.py
Print a list the chromosomes that are in a database (and their lengths)

#### combine_tracks.py
Combine counts from 2 or more tracks into a single new track. 
Very useful for combining data from multiple individuals, multiple replicates or from the 
forward and reverse strands.

#### liftover_track.py
Copy data in a track from one assembly to another (e.g. hg19 to hg18) 
using information from a UCSC liftover chain file. Here is an example of how data tracks 
can be copied from hg19 to hg18:

    # copy forward- and reverse-stranded data from hg19 to hg18
    python liftover_track.py hg19 hg18 ~/data/ucsc/hg19/liftover/hg19ToHg18.over.chain.gz \
        --rev_track uwdnase/wgEncodeUwDnaseWi38AlnRep2_rev uwdnase/wgEncodeUwDnaseWi38AlnRep2_fwd

    # copy unstranded data from hg19 to hg18
    python liftover_track.py hg19 hg18 /data/ucsc/hg19/liftover/hg19ToHg18.over.chain.gz uwdnase/wgEncodeUwDnaseWi38Aln

#### set_track_stats.py
Computes statistics for a track (n, mean, max, min, etc.) and stores them 
as attributes for each chromosome node in the HDF5 file. Attributes stored this way can be rapidly 
retrieved. You need to have write permissions for the track to run this script.

#### get_track_stats.py
Retrieves statistics for a track that have been pre-computed using set_track_stats.py


