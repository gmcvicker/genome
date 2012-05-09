
class TxtFileError(Exception):
    """An exception indicating a problem with a data table"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return str(self.value)


def read_rows(iterable, delim="\t"):
    """Creates an iterator, which can be used to read rows with named
    columns from a simple text file table"""
    iterator = iter(iterable)

    line_num = 1
    header = iterator.next().strip()
    while header.startswith("#"):
        line_num += 1
        header = iterator.next().strip()

    # parse the header line, get column names
    col_names = header.split(delim)
    n_col = len(col_names)
    idx2name = [None] * n_col
    idx = 0
    for col_name in col_names:
        idx2name[idx] = col_name
        idx += 1

    row_dict = {}
    
    # get the first line of data
    for line in iterator:
        line_num += 1
        line = line.strip()
        if not line.startswith("#"):
            cols = line.split(delim)
            if len(cols) != n_col:
                raise TxtFileError("Expected %d columns, but got %d "
                                   "on line %d" % (n_col, len(cols),
                                                   line_num))

            idx = 0
            for col in cols:
                row_dict[idx2name[idx]] = col
                idx += 1

            yield row_dict

    
