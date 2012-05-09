
import sys

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
        """returns a string representatin of this object"""
        return self.name

    def __cmp__(self, other):
        return cmp(self.idnum, other.idnum)


