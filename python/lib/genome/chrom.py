import gzip
import sys
import re

class Chromosome(object):
    def __init__(self, idnum=None, name=None, length=None,
                 is_sex=False, is_rand=False, is_hap=False,
                 is_mito=False, is_x=False, is_y=False, is_auto=False):
        self.idnum = idnum
        self.name = name
        self.length = length

        # is this a "random" chromosome?
        self.is_rand = is_rand
        
        # is this a sex chromosome?
        self.is_sex = is_sex

        # is this an alternate haplotype chromosome?
        self.is_hap = is_hap

        # is this the mitochondrial chromosome?
        self.is_mito = is_mito

        # is this the X chromosome
        self.is_x = is_x

        # is this the Y chromosome
        self.is_y = is_y

        # is this an autosome
        self.is_auto = is_auto
 
    def copy(self):
        """Creates a new chromosome object with the same attributes
        as this one"""
        return Chromosome(idnum=self.idnum,
                          name=None, length=None,
                          is_sex=False, is_rand=False, is_hap=False,
                          is_mito=False, is_x=False, is_y=False,
                          is_auto=False)
                          
    def __str__(self):
        """returns a string representation of this object"""
        return self.name

    def __cmp__(self, other):
        return cmp(self.idnum, other.idnum)








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
        
        chrom = Chromosome(name=words[0], length=words[1])
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




if __name__ == "__main__":

    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s chromInfo.txt.gz\n" % sys.argv[0])
        exit(2)
        
    chrom_list = parse_chromosomes(sys.argv[1])

    for chrom in chrom_list:
        sys.stderr.write(str(chrom) + "\n")
