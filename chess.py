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

def base_cur_r(x):
	return (base_top_r - r)/(base_h * base_h) * x * x + r

args = get_args()

origin = (0,0,0)
o_v = Vector(origin)

rq = 0.5 # ratio pawn base top to bottom radius
h = float(args.height)
r = float(args.radius)
th = float(args.thickness)
n = int(args.num_parts)
vn = int(args.num_vparts)
piece = args.piece
print("h=" + str(h))

def mk_file_path(orig_fname, suffix):
	parts = orig_fname.split(".")
	if len(parts) < 2:
		return orig_fname + "-" + suffix + ".stl"
	return parts[0] + "-" + suffix + "." + parts[1]

def make_pawn_bottom():
	x = 0.0
	dx = base_h / float(vn)
	verts = []
	faces = []
	cur_r = r

	for v_ind in range(0,vn):
		cur_r = base_cur_r(x)
		# print("x=" + str(x) + ", cur_r = " + str(cur_r))
		cur_verts = mu.get_circle_verts((0,0,x), cur_r, n)
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

	print("x="+str(x) + " dx=" + str(dx))
	#print("n_verts="+str(len(verts)))
	#print(faces)
	faces.append(range(0,n))
	faces.append(range((vn - 1) * n, vn * n))
	mu.create_mesh_from_data('Pawn', o_v, verts, faces, True)
	bpy.ops.export_mesh.stl(filepath=mk_file_path(args.save,"bottom"),ascii=False)

base_top_r = r * rq
base_q = 0.2
base_h = h * base_q

make_pawn_bottom()
