"""
This module contains a class, CoordLifter that is used to convert coordinates
from one assembly to another. The class requires a liftover HDF5 table to be present. This
table can be created using the script make_liftover_tab.py, which is part
of the genome repo.
"""

import sys
import genome.db
from genome.coord import CoordError


class CoordLifter(object):
    """A coordinate lifter object is used to convert coordinates from
    one assembly to another"""
    
    def __init__(self, from_assembly, to_assembly):        
        self.from_assembly = from_assembly
        self.to_assembly = to_assembly

        self.from_gdb = genome.db.GenomeDB(assembly=from_assembly)
        self.to_gdb = genome.db.GenomeDB(assembly=to_assembly)

        # lookup table of chromosomes in the "to" database, keyed on ID
        self.id_to_chrom = dict((x.idnum, x) for x in self.to_gdb.get_all_chromosomes())

        track_name = "liftover_%s_to_%s_tab" % (self.from_assembly, self.to_assembly)

        if not self.from_gdb.has_track(track_name):
            sys.stderr.write("liftover track named %s does not exist in genome db for "
                             "assembly %s\n" % (track_name, from_assembly))
        
        self.liftover_track = self.from_gdb.open_track(track_name)



    def convert(self, from_chrom, from_pos, from_strand=0):
        """Converts a position on the old assembly to a position on the new assembly.
        Returns a (chromosome object, position, strand) or None if the position could not
        be converted.
        """
        tab = self.liftover_track.h5f.getNode("/%s" % str(from_chrom))
        
        if from_pos > tab.nrows:
            raise CoordError("from_pos (%d) is greater than length of chromosome (%d) "
                             "in liftover table" % (from_pos, tab.nrows))

        row = tab[from_pos-1]
        to_chrom_id = row[0]
        to_chrom_pos = row[1]
        to_chrom_strand = row[2]
        
        if to_chrom_id == -1:
            return (None, None, None)
        
        to_chrom = self.id_to_chrom[to_chrom_id]
        
        return (to_chrom, to_chrom_pos, to_chrom_strand * from_strand)
    
        

        
        
if __name__ == "__main__":
    lifter = CoordLifter("hg18", "hg19")

    # chr1:995669 should convert to chr1:1005806
    new_coord = lifter.convert("chr1", 995669, 1)
    print(new_coord)

    # chr1:1007060 should convert to chr1:1017197
    new_coord = lifter.convert("chr1", 1007060, 1)
    print(new_coord)

    # chr22:15455353 should convert to chr22:17075353
    new_coord = lifter.convert("chr22", 15455353, 1)
    print(new_coord)
    
    # chr2:1203295 should be deleted in new
    new_coord = lifter.convert("chr2", 1203295, 1)    
    print(new_coord)
    
    

    
    
    
