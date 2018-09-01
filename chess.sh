#! /bin/bash

function die
{
  echo "$1"
  exit 1
}

OUT_FILE=${OUT_FILE:-pawn.stl}
PAWN_PARTS=3

# real pawn from chess set: radius = 8
# height = 30, top_r = 5
# base_h = 8
blender -b -P chess.py -- --save $OUT_FILE --height 300 --thickness 3 --num-parts 70 \
	--num-vparts 200 --radius 90 --edge-w 15 --piece pawn
