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
	parser.add_argument('-m', '--save', help="output file", required=True)
	parser.add_argument('-H', '--height', help="height", required=True)
	parser.add_argument('-f', '--flap-height', help="flap-height", required=True)
	parser.add_argument('-W', '--width', help="width", required=True)
	parser.add_argument('-w', '--flap-width', help="flap-width", required=True)
	parser.add_argument('-T', '--thickness', help="thickness", required=True)
	parser.add_argument('-t', '--flap-thickness', help="flap-thickness", required=True)
	parser.add_argument('-g', '--gap', help="gap", required=True)
	parsed_script_args, _ = parser.parse_known_args(script_args)
	return parsed_script_args


args = get_args()
h =  float(args.height)
flap_h = float(args.flap_height)
w =  float(args.width)
flap_w = float(args.flap_width)
th = float(args.thickness)
flap_th = float(args.flap_thickness)
gap = float(args.gap)

origin = (0,0,0)
mu.remove_cube()

inner_brick_verts = mu.get_brick_verts(origin, h, w, th)
inner_brick_faces = mu.get_hollow_brick_faces()
outer_brick_verts = mu.get_brick_verts(origin, h, w + gap, th + gap)
outer_brick_faces = mu.get_hollow_brick_faces()

mu.create_mesh_from_data('Inner', origin, inner_brick_verts, inner_brick_faces, False)
mu.create_mesh_from_data('Outer', origin, outer_brick_verts, outer_brick_faces, False)
bpy.ops.export_mesh.stl(filepath=args.save, ascii=False)


