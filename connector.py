import argparse
import bpy
import mathutils
import math
from mathutils import Vector
import sys
from math import radians

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
	parser.add_argument('-a', '--flap-angle', help='flap angle', required=True)
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
flap_angle = radians(float(args.flap_angle))
flap_delta_h = flap_th/math.tan(flap_angle)

origin = Vector((0,0,0))
mu.remove_cube()

inner_brick_verts = mu.get_brick_verts(origin, h, w, th)
inner_brick_faces = [list(reversed(f)) for f in mu.get_hollow_brick_faces()]
outer_brick_verts = mu.get_brick_verts(origin, h , w + gap, th + gap)
outer_brick_faces = mu.get_hollow_brick_faces()
gap_verts = inner_brick_verts + outer_brick_verts
gap_faces = []
next_gap_h = {5:7, 7:6, 6:4, 4:5, 3:1, 1:0, 0:2, 2:3}

for i in range(0,8):
	next_inner_i = next_gap_h[i]
	gap_faces.append([i, next_inner_i, next_inner_i + 8, i + 8])

# now some formulas resuling from solving the geometry
flap_h_long = flap_h + flap_delta_h
flap_h_part = flap_h_long - flap_h / 2.0
d_center = math.sqrt((flap_th * flap_th)/4.0 + flap_h_part * flap_h_part)
# angle between the diagonal and the lower side

d_angle_low = math.atan2(flap_th/2.0, flap_h/2.0 + flap_delta_h)
print("d_angle_low: " + str(d_angle_low))
lift_h = d_center * math.sin(flap_angle - d_angle_low);
print("lift_h: " + str(lift_h))
print("d_trap: " + str(d_center))
print("flap_delta_h: " + str(flap_delta_h))
flap_origin = Vector((0, 0, th/2 + lift_h))
#print(gap_verts)
#print(gap_faces)
#print(inner_brick_faces)
print(flap_angle)
flap_verts = mu.get_tp_prism_verts(origin, flap_h, flap_w, flap_th, flap_delta_h,
																	 0, math.pi - flap_angle, 0)
print(flap_verts)
flap_faces = mu.get_brick_faces()
mu.create_mesh_from_data('Inner', origin, inner_brick_verts, inner_brick_faces, False)
mu.create_mesh_from_data('Outer', origin, outer_brick_verts, outer_brick_faces, False)
mu.create_mesh_from_data('Gap sides', origin, gap_verts, gap_faces, False)
mu.create_mesh_from_data('Flap', flap_origin, flap_verts, flap_faces, False)
bpy.ops.export_mesh.stl(filepath=args.save, ascii=False)


