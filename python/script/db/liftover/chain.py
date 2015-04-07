
import gzip

import genome.coord
from genome.coord import Coord, CoordError
import genome.db



class Chain(object):
    def __init__(self, line, from_chrom_dict, to_chrom_dict):
        words = line.rstrip().split()

        #
        # Parse from coordinates for the this chain, being careful
        # that they look correct
        #
        from_chrom_name = words[2]
        from_chrom_size = int(words[3])

        if from_chrom_name in from_chrom_dict:
            from_chrom = from_chrom_dict[from_chrom_name]
        else:
            raise CoordError("chromosome '%s' does not exist "
                             "in 'from' database'" % from_chrom_name)

        if from_chrom.length != from_chrom_size:
            raise CoordError("chromosome length mismatch between chain "
                             "file and 'from' database %s:1-%d != %s:1-%d"
                             % (from_chrom_name, from_chrom_size,
                                from_chrom_name, from_chrom.length))

        if words[4] == '+':
            from_strand = 1
        else:
            raise CoordError("expected from strand from be '+'")
        
        from_start = int(words[5]) + 1
        from_end   = int(words[6])
        self.from_coord = Coord(from_chrom, from_start, from_end,
                                strand=from_strand)

        #
        # now parse "to" coordaintes
        #
        to_chrom_name = words[7]                
        to_chrom_size = int(words[8])

        if to_chrom_name in to_chrom_dict:
            to_chrom = to_chrom_dict[to_chrom_name]
        else:
            raise CoordError("chromosome %s does not exist "
                             "in 'to' database'" % to_chrom_name)

        if to_chrom.length != to_chrom_size:
            raise CoordError("chromosome length mismatch between chain "
                             "file and 'to' database %s:1-%d != %s1-%d"
                             % (to_chrom_name, to_chrom_size,
                                to_chrom_name, to_chrom.length))

        if words[9] == '+':
            to_strand = 1
        elif words[9] == '-':
            to_strand = -1
        else:
            raise CoordError("expected 'to' strand to be '+' or '-', "
                             "not '%s'" % words[9])

        if to_strand == 1:
            to_start = int(words[10]) + 1
            to_end = int(words[11])
        else:
            # coordinates in file are on reverse strand
            # and are 0-based half-open
            to_start = to_chrom_size - int(words[11]) + 1
            to_end = to_chrom_size - int(words[10])

        self.to_coord = Coord(to_chrom, to_start,
                              to_end, strand=to_strand)

        self.ori = to_strand

        self.blocks = []


    def add_block(self, block):
        # sanity check:  start / end should fall within chain bounds
        if (block.from_start < self.from_coord.start) or \
            (block.from_end > self.from_coord.end):
            raise CoordError("expected from block (%d-%d) to be within chain "
                             "from coords: %d-%d" % (block.from_start, block.from_end,
                                                     self.from_start, self.from_end))
        if (block.to_start < self.to_coord.start) or \
            (block.to_end > self.to_coord.end):
            raise CoordError("expected to block (%d-%d) to be within chain "
                             "from coords: %d-%d" % (block.to_start, block.to_end,
                                                     self.to_start, self.to_end))

        self.blocks.append(block)

    def __str__(self):
        return "%s %d %d %d %s %d %d" % (self.from_coord.chrom.name,
                                         self.from_coord.start,
                                         self.from_coord.end,
                                         self.ori,
                                         self.to_coord.chrom.name,
                                         self.to_coord.start,
                                         self.to_coord.end)


class ChainBlock(object):
    def __init__(self, from_start, from_end, to_start, to_end):
        self.from_start = from_start
        self.from_end = from_end
        self.to_start = to_start
        self.to_end = to_end

    def __str__(self):
        return "%d %d => %d %d" % (self.from_start, self.from_end,
                                   self.to_start, self.to_end)
        

        

def read_chain_file(chain_file, from_gdb, to_gdb):
    # read chromosomes from both databases
    from_chrom_dict = from_gdb.get_chromosome_dict()
    to_chrom_dict = to_gdb.get_chromosome_dict()

    if chain_file.endswith(".gz"):
        f = gzip.open(chain_file)
    else:
        f = open(chain_file)

    chain_list = []
        
    for line in f:
        if line.startswith('chain'):
            # this is the start of a new chain
            chain = Chain(line, from_chrom_dict, to_chrom_dict)
            # set starts (or ends) of block to the edge of this chain
            from_start = chain.from_coord.start
            from_end = None
            if chain.ori == 1:
                to_start = chain.to_coord.start
                to_end = None
            else:
                to_end = chain.to_coord.end
                to_start = None

            chain_list.append(chain)
            
        else:
            # this is a line describing a block within a chain
            # update ends of blocks
            words = line.rstrip().split()

            if len(words) == 0:
                # there is a blank line between chain lines
                continue
            
            size = int(words[0])
            from_end = from_start + size - 1

            if chain.ori == 1:
                to_end = to_start + size - 1
            else:
                # on reverse strand: to coords go other direction
                to_start = to_end - size + 1

            # create the block and add it to the chain
            block = ChainBlock(from_start, from_end, to_start, to_end)            
            chain.add_block(block)
            
            if len(words) == 3:
                # there are gaps between blocks.
                # update the start of the _next_ block to be past the end
                # of this one by a defined offset
                from_offset = int(words[1])
                to_offset = int(words[2])

                from_start = from_end + from_offset + 1

                if chain.ori == 1:
                    to_start = to_end + to_offset + 1
                else:
                    to_end = to_start - to_offset - 1
                    
            elif len(words) == 1:
                # this was the last line in chain, contains single ungapped alignment block
                # which we have already created above
                
                # sanity check: we should be at the end of the chain
                if from_end != chain.from_coord.end:
                    raise CoordError("end of last block (%d) does not "
                                     "match end of chain: %d" %
                                     (from_end, chain.from_coord.end))
                
                if chain.ori == 1:
                    if to_end != chain.to_coord.end:
                        raise CoordError("end of last block (%d) does not "
                                         "match end of chain: %d" %
                                         (to_end, chain.to_coord.end))
                else:
                    if to_start != chain.to_coord.start:
                        raise CoordError("start of last block (%d) does not "
                                         "match start of chain: %d" %
                                         (to_start, chain.to_coord.start))
            else:
                raise ValueError("expected line to have 1 or 3 tokens")

    f.close()

    return chain_list




if __name__ == "__main__":
    # for testing, read a chain file and print out blocks
    gdb18 = genome.db.GenomeDB(assembly="hg18")
    gdb19 = genome.db.GenomeDB(assembly="hg19")

    liftover_path = "/home/gmcvicker/data/ucsc/hg18/liftover/hg18ToHg19.over.chain.gz"
    
    chain_list = read_chain_file(liftover_path, gdb18, gdb19)

    for c in chain_list:
        print(str(c))

        for blk in c.blocks:
            print("  " + str(blk))
            
    
        
    
