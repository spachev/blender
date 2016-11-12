import argparse
import bpy
import mathutils
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
	parser.add_argument('-m', '--save', help="output file")
	parsed_script_args, _ = parser.parse_known_args(script_args)
	return parsed_script_args


args = get_args()
origin = Vector((0,0,0))
[x,y,z] = [0.707107, 0.258819, 0.965926]
verts = [[x,x,-1], [x,-x,-1], [-x,-x,-1], [-x,x,-1], [0,0,1], [1,1,1]]
#faces = [[1,0,4], [4,2,1], [4,3,2], [4,0,3], [0,1,2,3], [1,4,5],
#	[5,1,0], [4,5,0]]
#faces=[[0,1,2,3]]
verts = mu.mul_v(verts, 20)
print(verts)

cone1 = mu.create_mesh_from_data('DataCone', origin, verts, faces)

# export scene
bpy.ops.export_mesh.stl(filepath=args.save)
