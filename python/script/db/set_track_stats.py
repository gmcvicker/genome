
import sys
import numpy as np

import genome.db

import genome.trackstat as trackstat

def main():
    gdb = genome.db.GenomeDB()

    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <track_name>\n" % sys.argv[0])
        exit(2)

    track_name = sys.argv[1]

    track = gdb.open_track(track_name, "a")
    track_stat = trackstat.set_stats(gdb, track)
    sys.stderr.write("combined %s\n" % str(track_stat))
    track.close()


if __name__ == "__main__":
    main()
