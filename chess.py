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

class Piece:
	def __init__(self):
		pass

	def get_fname(self):
		return self.get_piece_name().lower() + ".stl"

	def save_part(self, part_name):
		bpy.ops.export_mesh.stl(filepath=mk_file_path(self.get_fname(), part_name),ascii=False)

	def make_ring_part(self, start_x, part_h, vn_arg, cur_outer_r_f, cur_inner_r_f,
										 part_name, all_in_one = False):
		if not all_in_one:
			mu.reset_scene()
		x = start_x
		dx = part_h / float(vn_arg)
		verts = []
		faces = []
		print("start_x="+str(x) + " dx=" + str(dx))

		for v_ind in range(0, vn_arg + 1):
			cur_inner_r = cur_inner_r_f(x)
			cur_outer_r = cur_outer_r_f(x)
			#print("x=" + str(x) + ", cur_inner_r = " + str(cur_inner_r) +
			#	", cur_outer_r = " + str(cur_outer_r))
			cur_inner_verts = mu.get_circle_verts((0,0,x), cur_inner_r, n)
			cur_outer_verts = mu.get_circle_verts((0,0,x), cur_outer_r, n)
			#cur_n_verts = len(verts)
			verts += cur_outer_verts + cur_inner_verts
			x += dx

			if v_ind > 0:
				outer_base_low = (v_ind - 1) * 2 * n
				inner_base_low = outer_base_low + n
				for i in range(0, n):
					outer_next_i = outer_base_low + (i + 1) % n
					inner_next_i = inner_base_low + (i + 1) % n
					#print(base_low)
					outer_face = [ outer_base_low + i, outer_next_i, outer_next_i + 2 * n,
								outer_base_low + i + 2 * n ]
					faces.append(outer_face)
					faces.append([ inner_base_low + i + 2 * n, inner_next_i + 2 * n,
											 inner_next_i, inner_base_low + i])

		print("final x="+str(x) + " final r = " + str(cur_outer_r))
		print("n_verts="+str(len(verts)))
		#print(faces)
		if cur_outer_r:
			outer_base = vn_arg * 2 * n
			inner_base = outer_base + n
			for i in range(0, n):
				i_next = (i + 1) % n
				faces.append([outer_base + i, outer_base + i_next,
										 inner_base + i_next, inner_base + i ])
		mu.create_mesh_from_data(self.get_piece_name() + ' ' + part_name, o_v, verts, faces, True)
		if not all_in_one:
			self.save_part(part_name)

	def make_conic_part(self, start_x, part_h, vn_arg, cur_r_f, part_name, all_in_one = False):
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
		mu.create_mesh_from_data(self.get_piece_name() + ' ' + part_name, o_v, verts, faces, True)
		if not all_in_one:
			self.save_part(part_name)

# three-part part piece
class Piece3(Piece):
	def __init__(self):
		super().__init__()
		self.base_h = h * self.base_q()
		self.middle_h = h * self.middle_q()
		self.top_h = h * self.total_q() - self.base_h - self.middle_h
		self.rq = 0.5 # ratio piece base top to bottom radius
		self.base_top_r = r * self.rq
		self.bend_h = self.base_h + self.middle_h * self.bend_q()
		self.bend_r = self.base_top_r * self.bend_rq()
		self.vn_base = int(vn * self.base_q())
		self.vn_middle = int(vn * self.middle_q())
		self.vn_top = int(vn * self.top_q())
		self.middle_top_r = self.middle_cur_r(self.middle_h + self.base_h)

	def bend_rq(self):
		return 0.6
	def bend_q(self):
		return 0.7
	# used for stretching piece
	# in particular pawns should be shorter than other pieces
	def total_q(self):
		return 1.0
	def top_q(self):
		return self.total_q() - self.base_q() - self.middle_q()
	def base_cur_r(self, x):
		return (self.base_top_r - r)/(self.base_h * self.base_h) * x * x + r
	def middle_cur_r(self, x):
		if x > self.bend_h:
			return self.bend_r
		return (self.base_top_r - self.bend_r) / \
			((self.base_h - self.bend_h) * \
			(self.base_h - self.bend_h)) * (x - self.bend_h) \
		* (x - self.bend_h) + self.bend_r

	def make_base(self, all_in_one):
		self.make_conic_part(0, self.base_h, self.vn_base, self.base_cur_r,
												"base", all_in_one)

	def make_top(self, all_in_one):
		self.make_conic_part(self.middle_h + self.base_h, self.top_h, self.vn_top,
												self.top_cur_r, "top", all_in_one)

	def make_middle(self, all_in_one):
		self.make_conic_part(self.base_h, self.middle_h, self.vn_middle,
												self.middle_cur_r, "middle", all_in_one)

	def make(self):
		for all_in_one in [False, True]:
			mu.reset_scene()
			self.make_base(all_in_one)
			self.make_middle(all_in_one)
			self.make_top(all_in_one)
		bpy.ops.export_mesh.stl(filepath=self.get_fname(), ascii=False)

# real pawn: h = 30
class Pawn(Piece3):
	def __init__(self):
		super().__init__()

		# top_shift_h + top_r = top_h
		# top_shift_h ^ 2 + middle_top_r ^ 2 = top_r ^ 2
		# (top_h - top_r) ^ 2 + middle_top_r ^ 2 = top_r ^ 2
		# (top_h - top_r) ^ 2 = top_r ^ 2 - middle_top_r ^ 2
		# top_h ^ 2 - 2 * top_h * top_r + top_r ^2 = top_r ^2 - middle_top_r ^ 2
		# top_h ^ 2 - 2 * top_h * top_r = - middle_top_r ^ 2
		# top_r = ( top_h ^ 2  + middle_top_r ^ 2)/ (2*top_h)

		self.top_r = ( self.top_h * self.top_h  + \
			self.middle_top_r * self.middle_top_r) / (2.0 * self.top_h)
		self.top_shift_h = self.top_h - self.top_r

		print("base_h=" + str(self.base_h) + " bend_h = " + str(self.bend_h))
		print("bend_r=" + str(self.bend_r))
		print("top_shift_h = " + str(self.top_shift_h) + " top_h = " + str(self.top_h))

	def get_piece_name(self):
		return "Pawn"
	def base_q(self):
		return 0.25
	def middle_q(self):
		return 0.4

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

# real bishop: top_max_r = 6
# bend_r = 4
# h = 48
class Bishop(Piece3):
	def __init__(self):
		super().__init__()
		self.sin_q = 2.0
		# joint point of sinusoides
		self.sin_joint_h =  self.top_h *  1.0/(self.sin_q + 1.0)
		print("sin_joint_h="+str(self.sin_joint_h))
		self.top_max_rq = 1.5
		self.top_max_r = self.bend_r * self.top_max_rq
		self.sin_k = math.pi / (2.0*self.sin_joint_h)
		self.shift_top_x = math.asin(self.bend_r / self.top_max_r) / self.sin_k
		self.top_h = h * (self.total_q() - self.base_q() - \
			self.middle_q()) - self.shift_top_x

	def top_cur_r(self,x):
		x -= self.base_h + self.middle_h - self.shift_top_x
		#print("x=" + str(x))
		if x < self.sin_joint_h:
			return self.top_max_r * math.sin(self.sin_k * x)
		return self.top_max_r * math.cos(math.pi*(x - \
			self.sin_joint_h)/(2.0 * self.sin_joint_h * self.sin_q))

	def get_piece_name(self):
		return "Bishop"
	def total_q(self):
		return 1.2
	def base_q(self):
		return 0.3
	def middle_q(self):
		return 0.5

class RoyalPiece(Piece3):
	def __init__(self):
		super().__init__()
		self.funnel_h = self.top_h * self.funnel_q()
		self.funnel_r = self.bend_r * self.funnel_rq()
		self.crown_r = self.bend_r * self.crown_rq()
	def total_q(self):
		return 1.3
	def base_q(self):
		return 0.3
	def middle_q(self):
		return 0.7
	def funnel_q(self):
		return 0.5
	def funnel_rq(self):
		return 1.6
	def crown_rq(self):
		return 0.8
	def top_cur_r(self, x):
		x -= self.base_h + self.middle_h
		if x < self.funnel_h:
			return self.bend_r + x * (self.funnel_r - self.bend_r) / self.funnel_h
		return self.crown_cur_r(x)

class King(RoyalPiece):
	def __init__(self):
		super().__init__()
	def get_piece_name(self):
		return "King"
	def crown_cur_r(self, x):
		return self.crown_r * (self.top_h - x)/(self.top_h - self.funnel_h)

class Queen(RoyalPiece):
	def __init__(self):
		super().__init__()
		self.crown_r = self.top_h - self.funnel_h
	def get_piece_name(self):
		return "Queen"
	def funnel_q(self):
		return 0.7
	def crown_cur_r(self, x):
		x -= self.funnel_h
		return math.sqrt(self.crown_r * self.crown_r - x * x)

# real rook: h = 36
class Rook(Piece3):
	def __init__(self):
		super().__init__()
		self.top_r = self.bend_r * self.top_rq()
		self.pit_r = self.bend_r * (self.top_rq() - self.top_wall_rq())
		self.pit_h = self.top_h * self.top_pit_q()
		self.top_d_pit = self.top_h - self.pit_h
	def get_piece_name(self):
		return "Rook"
	def total_q(self):
		return 1.1
	def base_q(self):
		return 0.3
	def middle_q(self):
		return 0.6
	def top_rq(self):
		return 1.5
	def top_wall_rq(self):
		return 0.3
	def top_pit_q(self):
		return 0.4
	def bend_rq(self):
		return 0.8
	def top_cur_r(self, x):
		return self.top_r
	def pit_cur_r(self, x):
		return self.pit_r
	def make_top(self, all_in_one):
		if not all_in_one:
			mu.reset_scene()
		self.make_conic_part(self.base_h + self.middle_h, self.top_d_pit, self.vn_top,
			self.top_cur_r, "pitbottom", True)
		self.make_ring_part(self.base_h + self.middle_h + self.top_d_pit, self.pit_h,
			self.vn_top,
			self.top_cur_r, self.pit_cur_r, "pit", True)
		if not all_in_one:
			self.save_part("top")

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
bishop = Bishop()
bishop.make()
king = King()
king.make()
queen = Queen()
queen.make()
rook = Rook()
rook.make()
