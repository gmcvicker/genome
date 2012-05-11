#!/usr/bin/python

import genome.db

gdb = genome.db.GenomeDB()

for trackname in gdb.list_tracks():
    print trackname
