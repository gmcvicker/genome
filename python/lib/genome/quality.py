#!/usr/bin/env python


def qual_str_to_codes(qual_str):
    """Converts a string of space-delimited quality values to an
    ascii string representation, where the ascii value of each character
    is the quality value"""
    return "".join([chr(int(x)) for x in qual_str.split(" ")])


def qual_codes_to_str(qual_codes):
    """Converts ascii quality scores to a human-readible string of
    space-delimited numbers"""
    return " ".join([str(ord(x)) for x in qual_codes])



if __name__ == "__main__":
    qual_str = "45 50 60"
    ascii_str = qual_str_to_codes(qual_str)
    new_qual_str = qual_codes_to_str(ascii_str)
    
    print("original qual str: " + qual_str)
    print("ascii quals str: " + ascii_str)
    print("new qual str: " + new_qual_str)
