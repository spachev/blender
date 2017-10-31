#! /bin/bash

function die
{
  echo "$1"
  exit 1
}

OUT_FILE=${OUT_FILE:-pawn.stl}
PAWN_PARTS=3

blender -b -P chess.py -- --save $OUT_FILE --height 500 --thickness 3 --num-parts 70 \
	--num-vparts 20 --radius 120 --piece pawn
