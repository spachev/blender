#! /usr/bin/python

import argparse
import bpy
import mathutils
import math
from mathutils import Vector
import sys

sys.path.append('.')

import meshutil as mu

def get_args():
	parser = argparse.ArgumentParser()

	# get all script args
	_, all_arguments = parser.parse_known_args()
	double_dash_index = all_arguments.index('--')
	script_args = all_arguments[double_dash_index + 1: ]

	# add parser rules
	parser.add_argument('--save', help="output file", required=True)
	parser.add_argument('--width', help="width",  required=True)
	parser.add_argument('--height', help="height", required=True)
	parser.add_argument('--length', help="length", required=True)
	parser.add_argument('--num-parts', help="number of circle segements", required=True)
	parser.add_argument('--num-vparts', help="number of vertical segements", required=True)
	parsed_script_args, _ = parser.parse_known_args(script_args)
	return parsed_script_args

args = get_args()


origin = (0,0,0)
o_v = Vector(origin)

h = float(args.height)
w = float(args.width)
l = float(args.length)
