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
	parser.add_argument('-m', '--save', help="output file", required=True)
	parser.add_argument('-l', '--length', help="length",  required=True)
	parser.add_argument('-t', '--thickness', help="thickness", required=True)
	parser.add_argument('-n', '--num-parts', help="number of circle segements", required=True)
	parser.add_argument('-d', '--ball-diameter', help="ball diameter",  required=True)
	parser.add_argument('-D', '--ring-diameter', help="ring diameter",  required=True)
	parser.add_argument('-w', '--bar-width', help="bar-width",  required=True)
	parsed_script_args, _ = parser.parse_known_args(script_args)
	return parsed_script_args

def add_cylynder(name,o,h,r,n,rev_faces):
	top_center=(o[0],o[1],o[2]+h/2.0)
	bottom_center=(o[0],o[1],o[2]-h/2.0)
	top_circle_verts = mu.get_circle_verts(top_center,r,n)
	#top_circle_faces=[[i for i in range(0,n)]]
	bottom_circle_verts = mu.get_circle_verts(bottom_center,r,n)
	verts = top_circle_verts + bottom_circle_verts
	bottom_circle_faces=[[i+n for i in range(0,n)]]
	faces = [[i,i+1,i+1+n,i+n] for i in range(0,n-1)]
	faces = faces + [[n-1,0,n,2*n-1]]

	if rev_faces:
		for f in faces:
			f = f.reverse()

	mu.create_mesh_from_data(name, origin, verts, faces, False)
	return verts

def add_ring(name,o,h,r_in,r_out,n):
	in_verts = add_cylynder(name + " inner cylinder", o, h, r_in, n, False)
	out_verts = add_cylynder(name + " outer cylinder", o, h, r_out, n, True)
	verts = in_verts + out_verts
	top_faces = [[i+2*n,i+2*n+1,i+1,i] for i in range(0,n-1)]
	top_faces = top_faces + [[0,n-1,3*n-1,2*n]]
	bottom_faces = [[i+n,i+1+n,i+1+3*n,i+3*n] for i in range(0,n-1)]
	bottom_faces = bottom_faces + [[n,2*n-1,4*n-1,3*n]]
	mu.create_mesh_from_data(name + ' top face', origin, verts,
													 top_faces + bottom_faces, False)

args = get_args()
origin = (0,0,0)
th = float(args.thickness)
l = float(args.length)
n = int(args.num_parts)
ball_d = float(args.ball_diameter)
ring_d = float(args.ring_diameter)
w = float(args.bar_width)

bar_l = l - 2.0 * ring_d
brick_verts = mu.get_brick_verts(origin,bar_l + w/2.0,w,th)
brick_faces = mu.get_brick_faces()
mu.remove_cube()
mu.create_mesh_from_data('Bar', origin, brick_verts, brick_faces, False)
left_ring_center=(bar_l/2.0+ring_d/2.0,0,0)
right_ring_center=(-(bar_l/2.0+ring_d/2.0),0,0)
add_ring('left ring', left_ring_center, th, ball_d/2.0, ring_d/2.0, n)
add_ring('right ring', right_ring_center, th, ball_d/2.0, ring_d/2.0, n)

#print(brick_verts)
bpy.ops.export_mesh.stl(filepath=args.save,ascii=False)




