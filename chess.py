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


def mk_file_path(orig_fname, suffix):
	parts = orig_fname.split(".")
	if len(parts) < 2:
		return orig_fname + "-" + suffix + ".stl"
	return parts[0] + "-" + suffix + "." + parts[1]

class Pawn:
	def __init__(self):
		self.rq = 0.5 # ratio pawn base top to bottom radius
		self.bend_rq = 0.6 # ratio of bend_top_r to base_top_r
		self.base_top_r = r * self.rq
		self.base_q = 0.25
		self.middle_q = 0.4
		self.bend_q = 0.7
		self.top_q = 1.0 - (self.base_q + self.middle_q)
		self.base_h = h * self.base_q
		self.middle_h = h * self.middle_q
		self.bend_h = self.base_h + self.middle_h * self.bend_q
		self.bend_r = self.base_top_r * self.bend_rq
		self.middle_top_r = self.middle_cur_r(self.middle_h + self.base_h)

		# top_shift_h + top_r = top_h
		# top_shift_h ^ 2 + middle_top_r ^ 2 = top_r ^ 2
		# (top_h - top_r) ^ 2 + middle_top_r ^ 2 = top_r ^ 2
		# (top_h - top_r) ^ 2 = top_r ^ 2 - middle_top_r ^ 2
		# top_h ^ 2 - 2 * top_h * top_r + top_r ^2 = top_r ^2 - middle_top_r ^ 2
		# top_h ^ 2 - 2 * top_h * top_r = - middle_top_r ^ 2
		# top_r = ( top_h ^ 2  + middle_top_r ^ 2)/ (2*top_h)

		self.top_h = h - self.base_h - self.middle_h
		self.top_r = ( self.top_h * self.top_h  + \
			self.middle_top_r * self.middle_top_r) / (2.0 * self.top_h)
		self.top_shift_h = self.top_h - self.top_r

		print("base_h=" + str(self.base_h) + " bend_h = " + str(self.bend_h))
		print("bend_r=" + str(self.bend_r))
		print("top_shift_h = " + str(self.top_shift_h) + " top_h = " + str(self.top_h))

		self.vn_base = int(vn * self.base_q)
		self.vn_middle = int(vn * self.middle_q)
		self.vn_top = int(vn * self.top_q)

	def base_cur_r(self, x):
		return (self.base_top_r - r)/(self.base_h * self.base_h) * x * x + r

	# f(base_h) = base_top_r
	# f(bend_h) = bend_r
	# f(x) = a * (x - bend_h)^2 + bend_r
	# a * (base_h - bend_h) ^ 2 + bend_r = base_top_r
	# a * (base_h - bend_h) ^ 2 = base_top_r - bend_r
	# a = (base_top_r - bend_r) / (base_h - bend_h) ^ 2
	def middle_cur_r(self, x):
		if x > self.bend_h:
			return self.bend_r
		return (self.base_top_r - self.bend_r) / \
			((self.base_h - self.bend_h) * \
			(self.base_h - self.bend_h)) * (x - self.bend_h) \
		* (x - self.bend_h) + self.bend_r

	# top_shift_h + top_r = top_h
	# (top_shift_h + top_r) = 2 * top_r * top_ball_q
	# top_h = 2 * top_r * top_ball_q
	# top_r = top_h /( 2 * top_ball_q)
	def top_cur_r(self, x):
		dx = x - self.top_shift_h - self.middle_h - self.base_h
		y_sq = self.top_r * self.top_r - dx * dx
		#print("y_sq = " + str(y_sq) + " x = " + str(x))
		if y_sq > 0.0:
			return math.sqrt(y_sq)
		else:
			return 0

	def make_pawn_part(self, start_x, part_h, vn_arg, cur_r_f, part_name, all_in_one = False):
		if not all_in_one:
			mu.reset_scene()
		x = start_x
		dx = part_h / float(vn_arg)
		verts = []
		faces = []
		print("start_x="+str(x) + " dx=" + str(dx))

		for v_ind in range(0, vn_arg + 1):
			cur_r = cur_r_f(x)
			#print("x=" + str(x) + ", cur_r = " + str(cur_r))
			cur_verts = mu.get_circle_verts((0,0,x), cur_r, n)
			cur_n_verts = len(verts)
			verts += cur_verts
			x += dx

			if v_ind > 0:
				base_low = (v_ind - 1) * n
				for i in range(0, n):
					next_i = base_low + (i + 1) % n
					#print(base_low)
					faces.append([ base_low + i, next_i, next_i + n,
								base_low + i + n ])

		print("final x="+str(x) + " final r = " + str(cur_r))
		print("n_verts="+str(len(verts)))
		#print(faces)
		faces.append(range(0,n))
		if cur_r:
			faces.append(range((vn_arg) * n, (vn_arg + 1) * n))
		mu.create_mesh_from_data('Pawn ' + part_name, o_v, verts, faces, True)
		if not all_in_one:
			bpy.ops.export_mesh.stl(filepath=mk_file_path(args.save, part_name),ascii=False)

	def make_pawn_base(self, all_in_one):
		self.make_pawn_part(0, self.base_h, self.vn_base, self.base_cur_r,
												"base", all_in_one)

	def make_pawn_top(self, all_in_one):
		self.make_pawn_part(self.middle_h + self.base_h, self.top_h, self.vn_top,
												self.top_cur_r, "top", all_in_one)

	def make_pawn_middle(self, all_in_one):
		self.make_pawn_part(self.base_h, self.middle_h, self.vn_middle,
												self.middle_cur_r, "middle", all_in_one)

	def make(self):
		for all_in_one in [False, True]:
			mu.reset_scene()
			self.make_pawn_base(all_in_one)
			self.make_pawn_middle(all_in_one)
			self.make_pawn_top(all_in_one)
		bpy.ops.export_mesh.stl(filepath='pawn.stl', ascii=False)

args = get_args()

origin = (0,0,0)
o_v = Vector(origin)

h = float(args.height)
r = float(args.radius)
th = float(args.thickness)
n = int(args.num_parts)
vn = int(args.num_vparts)
piece = args.piece

pawn = Pawn()
pawn.make()

