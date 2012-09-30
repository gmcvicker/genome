import subprocess
import os
import gzip


def check_open(filename, mode="r"):
    """Tries to open file and return filehandle. Takes into account
    that file may be gzipped. Raises exception if mode is write and 
    file already exists."""
    if mode.startswith("w") and os.path.exists(filename):
        raise IOError("file %s already exists" % filename)

    if filename.endswith(".gz"):
        if mode == "w":
            mode = "wb"
        elif mode == "r":
            mode = "rb"

        return gzip.open(filename, mode)

    return open(filename, mode)



def count_lines(filename):

    if not os.path.exists(filename):
        raise IOError("file '%s' does not exist" % filename)

    if not os.path.isfile(filename):
        raise IOError("'%s' is not a regular file" % filename)
    
    if filename.endswith('gz'):
        p1 = subprocess.Popen(['zcat', filename],
                              stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout,
                               stdout=subprocess.PIPE)
        
        p1.stdout.close()

        out = p2.communicate()[0]
    else:
        out = subprocess.Popen(['wc', '-l', filename],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT).communicate()[0]

    return int(out.split()[0])
