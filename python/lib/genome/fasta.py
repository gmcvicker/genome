#!/usr/bin/env python

def read_fasta(iterable, linejoin=""):
    iterator = iter(iterable)
    line = iterator.next()
    
    if not line.startswith(">"):
        raise ValueError("first line does not start with '>'")
    
    header = line[1:].rstrip()
    
    lines = []
    for line in iterator:
        if line.startswith(">"):
            # return from here, but restart here when we re-enter this function
            yield header, linejoin.join(lines)
            
            header = line[1:].rstrip()
            lines = []
        else:
            lines.append(line.rstrip())
    
    yield header, linejoin.join(lines)


def read_fasta_lines(iterable):
    iterator = iter(iterable)
    line = iterator.next()

    if not line.startswith(">"):
        raise ValueError("first line does not start with '>'")

    header = line[1:].rstrip()

    lines = []
    for line in iterator:
        if line.startswith(">"):
            # return from here, but restart here when we re-enter
            # this function
            yield header, lines

            lines = []
            header = line[1:].rstrip()
        else:
            lines.append(line.rstrip())

    yield header, lines


def write_fasta(file, id, seq_str, line_width=60):
    file.write(">" + id + "\n")
    for p in xrange(0, len(seq_str), line_width):
        file.write(seq_str[p:p+line_width] + "\n")



def write_vals(file, id, vals, line_width=60):
    file.write(">" + id + "\n")
    for p in xrange(0, len(vals), line_width):
        file.write(" ".join([str(q) for q in vals[p:p+line_width]]) + "\n")
