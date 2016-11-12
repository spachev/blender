#! /bin/bash
OUT_FILE=${OUT_FILE:-helicopter-part.stl}
blender -b -P helicopter_part.py -- -m $OUT_FILE -l 11.43 -D 4.19 -t 1.27 -w 1.27 -d 2.00 -n 70
