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
	parser.add_argument('-r', '--radius', help="radius",  required=True)
	parser.add_argument('-H', '--height', help="height", required=True)
	parser.add_argument('-k', '--radius-k', help="radius k", required=True)
	parser.add_argument('-t', '--thickness', help="thickness", required=True)
	parser.add_argument('-v', '--vthickness', help="vthickness", required=True)
	parser.add_argument('-n', '--num-parts', help="number of circle segements", required=True)
	parser.add_argument('-T', '--text', help="Text on the cup", required=True)
	parser.add_argument('-f', '--font-size', help="Font size", required=True)
	parsed_script_args, _ = parser.parse_known_args(script_args)
	return parsed_script_args

args = get_args()

origin = (0,0,0)
th = float(args.thickness)
h = float(args.height)
v_th = float(args.vthickness)
radius_k = float(args.radius_k)
font_size = int(args.font_size)
cup_text = args.text
out_bottom_center = (0,0,0)
bottom_center = (0,0,v_th)
top_center = (0,0,h)
out_top_center = (0,0,h)
o_v = Vector(origin)
bottom_r = float(args.radius)
n_parts = int(args.num_parts)
top_r = radius_k * bottom_r
bottom_verts = mu.get_circle_verts(bottom_center, bottom_r, n_parts)
top_verts = mu.get_circle_verts(top_center, top_r, n_parts)
out_bottom_verts = mu.get_circle_verts(out_bottom_center, bottom_r + v_th, n_parts)
ext_bottom_verts = mu.get_circle_verts(bottom_center, bottom_r + v_th, n_parts)
out_top_verts = mu.get_circle_verts(out_top_center, top_r + v_th, n_parts)

#print(verts)
bottom_face = [ i for i in range(0, n_parts)]
top_face = [ i + n_parts for i in range(0, n_parts)]
out_bottom_face = [ i + 2 * n_parts for i in range(n_parts-1,0,-1)]
out_top_face = [ i + 3 * n_parts for i in range(0, n_parts)]

faces = [bottom_face, out_bottom_face]
verts = bottom_verts + top_verts + out_bottom_verts + out_top_verts + \
	ext_bottom_verts

# inner wall faces
for i in range(0,n_parts-1):
	faces.append([i + 1, i, i+n_parts, i+n_parts+1])

faces.append([0, n_parts-1, n_parts-1+n_parts, n_parts])

# outer wall faces
for i in range(n_parts*2,n_parts*3-1):
	faces.append([i, i+1, i+1+n_parts, i+n_parts])

faces.append([n_parts*2, 3*n_parts-1, 3*n_parts-1+n_parts, 3*n_parts])

# connect outer top to top
for i in range(n_parts,n_parts*2-1):
	faces.append([i+1, i, i+2*n_parts, i+2*n_parts + 1])

faces.append([n_parts, 2*n_parts-1, 2*n_parts-1+2*n_parts, 3*n_parts])

# connect inner to outer bottom for re-enforcement
for i in range(0,n_parts-1):
	faces.append([i+1, i, i+4*n_parts, i+4*n_parts + 1])

faces.append([0, 4*n_parts-1, n_parts-1+4*n_parts, n_parts])

#print(faces)
#print(verts)

#me = mu.create_mesh_from_data('Cup', o_v, verts, faces, False)
text_r = (bottom_r + v_th)*(1.0 + radius_k)/2.0
text_ov = Vector((0, 0, h/2.0))
#new_o = Vector((0,0,2*h))
#verts = mu.get_circle_verts(new_o, bottom_r, n_parts)
#face1 = [ i for i in range(0,n_parts) ]
#faces = [face1]
mu.remove_cube()
mu.add_text(cup_text, text_ov, text_r, n_parts, font_size)
mu.create_mesh_from_data('Test', o_v, verts, faces, False)
#print(len(bpy.context.scene.objects))
#print(len(bpy.context.scene.objects))
bpy.ops.export_mesh.stl(filepath=args.save,ascii=False)
