import gzip
import sys
import re
import string


LINES_PER_RECORD = 4
MAX_LINE_LEN = 2048


def replace_non_printable(str):
    can_print = set(string.printable)

    new_chars = [x for x in str]
    
    for i in range(len(new_chars)):
        if new_chars[i] not in can_print:
            new_chars[i] = "\{" + str(ord(new_chars[i])) + "}"
    
    return "".join(new_chars)


    

class FastqError(Exception):
    def __init__(self, msg, lines=None):
        self.msg = msg
	self.lines = lines
        if lines:
          self.__fix_lines()
        

    def __fix_lines(self):
        for i in range(len(self.lines)):
            if len(self.lines[i]) > MAX_LINE_LEN:
                self.lines[i] = self.lines[i][0:MAX_LINE_LEN]
                self.lines[i] = replace_non_printable(self.lines[i])
 

    def __str__(self):
        if self.lines:
            return self.msg + ":\n  " + "\n  ".join(self.lines)
        else:
            return self.msg



def check_fastq(lines):
    if not lines[0].startswith("@"):
        raise FastqError("read header does not start with '@'", lines)

    if not lines[2].startswith("+"):
        raise FastqError("qual header does not start with '+'", lines)

    if len(lines[1]) != len(lines[3]):
        raise FastqError("qual len does not match read len", lines)

    if len(lines[2]) > 1:
        # if there is a complete header before the quals, verify it matches
        # the read header
        if lines[0][1:] != lines[2][1:]:
            raise FastqError("qual header does not match read header", lines)
    
        
                        
    
    
def read_fastq(filename):
    """Creates a fastq reader that returns lines at a time. Raises an
    exception if the format looks incorrect."""

    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)

    read_len = None
    line_num = 0
    lines = []

    for line in f:
        if line_num == LINES_PER_RECORD:
            check_fastq(lines)
            yield lines
            lines = []
            line_num = 0
        line_num += 1        
        lines.append(line.rstrip())

    if line_num != LINES_PER_RECORD:
        raise FastqError("number of lines in file was not multiple of 4", lines)

    check_fastq(lines)
    
    yield lines



def write_fastq(f, lines):
    # write lines to output files
    for i in range(LINES_PER_RECORD):
        f.write(lines[i] + "\n")
