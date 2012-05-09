
import sys
import tables

from genome.chrom import Chromosome

from adaptor import DBAdaptor


class ChromosomeAdaptor(DBAdaptor):
    def __init__(self, chrom_track):
        super(ChromosomeAdaptor, self).__init__(chrom_track)
    
    def fetch_by_name(self, name_container=[]):
        """Returns a chromosome iterator for the provided GenomeDB. Only
        chromosomes with names in the provided container are yielded by
        the iterator."""

        if type(name_container) is str:
            # wrap string in list because we
            # don't want statements like 'name in container'
            # to be True for "chr2" and "chr22", for example
            name_container = [name_container]

        chrom = Chromosome()

        for row in self.track.h5f.root.chromosome:
            if row['name'] in name_container:
                chrom.idnum = row['idnum']
                chrom.is_auto = row['is_auto']
                chrom.is_hap = row['is_hap']
                chrom.is_mito = row['is_mito']
                chrom.is_sex = row['is_sex']
                chrom.is_x = row['is_x']
                chrom.is_y = row['is_y']
                chrom.length = row['length']
                chrom.name = row['name']

                yield chrom


    def fetch_by_type(self, get_rand=False,
                      get_auto=True, get_sex=True, get_x=True,
                      get_y=False, get_hap=False, get_mito=False):
        """Returns a chromosome iterator for the provided GenomeDB."""

        chrom = None

        flags = {'is_rand' : get_rand,
                 'is_auto' : get_auto,
                 'is_sex' : get_sex,
                 'is_x' : get_x,
                 'is_y' : get_y,
                 'is_hap' : get_hap,
                 'is_mito' :  get_mito}

        conditions = []
        for k, v in flags.items():
            if not v:
                # don't get this chromosome type
                conditions.append("(%s == False)" % k)

        row_iter = None
        if len(conditions) > 0:
            query = " & ".join(conditions)
            row_iter = self.track.h5f.root.chromosome.where(query)
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

            yield chrom


    def fetch_name_dict(self):
        """Retrieves chromosomes from the provided GenomeDB,
        and returns a dictionary keyed on chromosome name"""

        chrom_dict = {}
        for row in self.track.h5f.root.chromosome:
            chrom = Chromosome(idnum=row['idnum'],
                               name=row['name'],
                               length=row['length'],
                               is_auto=row['is_auto'],
                               is_hap=row['is_hap'],
                               is_mito=row['is_mito'],
                               is_sex=row['is_sex'],
                               is_x=row['is_x'],
                               is_y=row['is_y'])

            if chrom.name in chrom_dict:
                raise ValueError("duplicate chromosome name: %s" % chrom.name)
            chrom_dict[chrom.name] = chrom

        return chrom_dict



if __name__ == "__main__":
    # read configuration file
    gdb = GenomeDB()
    
    chrom_track = self.open_track("chromosome")
    chrom_adp = ChromosomeAdaptor(chrom_track)
        
    for chrom in chrom_adp.fetch_by_type(db):
        print "%s ID:%d len:%d" % (chrom.name, chrom.idnum, chrom.length)
