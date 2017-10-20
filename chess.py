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
	parser.add_argument('--radius', help="radius",  required=True)
	parser.add_argument('--height', help="height", required=True)
	parser.add_argument('--thickness', help="thickness", required=True)
	parser.add_argument('--num-parts', help="number of circle segements", required=True)
	parser.add_argument('--num-vparts', help="number of vertical segements", required=True)
	parser.add_argument('--piece', help='chess piece name', required=True)
	parsed_script_args, _ = parser.parse_known_args(script_args)
	return parsed_script_args

args = get_args()

origin = (0,0,0)
o_v = Vector(origin)

rq = 0.5 # ratio pawn base height to the radius
h = float(args.height)
r = float(args.radius)
th = float(args.thickness)
n = int(args.num_parts)
vn = int(args.num_vparts)
piece = args.piece
print("h=" + str(h))
x = 0
dx = rq * r / float(vn)
verts = []
faces = []
cur_r = r

for v_ind in range(0,vn):
	cur_r = cur_r - x
	cur_verts = mu.get_circle_verts((0,0,cur_r - x), cur_r, n)
	cur_n_verts = len(verts)
	verts += cur_verts
	x += dx

	if v_ind > 0:
		for i in range(0, n):
			base_low = (v_ind - 1) * n
			next_i = base_low + (i + 1) % n
			#print(base_low)
			faces.append([ base_low + i, next_i, next_i + n,
						base_low + i + n ])

#print("n_verts="+str(len(verts)))
#print(faces)
mu.create_mesh_from_data('Pawn', o_v, verts, faces, True)
bpy.ops.export_mesh.stl(filepath=args.save,ascii=False)
