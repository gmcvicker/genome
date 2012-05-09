import subprocess
import os

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
