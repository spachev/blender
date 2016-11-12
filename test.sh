#! /bin/bash
OUT_FILE=${OUT_FILE:-cup1.stl}
blender -b -P cup.py -- -m $OUT_FILE -T William -H 100 -r 30 -k 1.3 -t 3 -v 4 -n 70 -f 20
