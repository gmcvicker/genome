import re

def from_file(filename, skip_header=False, flip_key_val=False):
    """Reads a file containing tab-separated key/value pairs and
    returns the resulting dictionary. Values in the returned dictionary
    are lists, that may contain multiple elements"""
    f = open(filename, "r")
    if skip_header:
        f.readline()

    d = {}    
    for line in f:
        words = line.strip().split()
        if len(words) == 2:
            if flip_key_val:
                k = words[1]
                v = words[0]
            else:
                k = words[0]
                v = words[1]

            if k in d:
                d[k].append(v)
            else:
                d[k] = [v]
        else:
            raise ValueError("Too many words on line:\n%s" % line)
            
    f.close()

    return d



def merge(dict_list):
    """Takes multiple dictionaries that have lists for values, and
    returns a new combined dictionary"""
    new_dict = {}

    for d in dict_list:
        for k, v in d.items():
            if k in new_dict:
                new_dict[k].append(v)
            else:
                new_dict[k] = v

    return new_dict
