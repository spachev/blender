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
	parser.add_argument('--edge-w', help="edge width", required=True)
	parser.add_argument('--num-parts', help="number of circle segements", required=True)
	parser.add_argument('--num-vparts', help="number of vertical segements", required=True)
	parser.add_argument('--piece', help='chess piece name', required=False)
	parser.add_argument('--part', help='chess piece part', default="all", required=False)
	parser.add_argument('--start-height', help='start height for partial print', default=0, required=False)
	parser.add_argument('--end-height', help='end height for partial print', default=0, required=False)
	parsed_script_args, _ = parser.parse_known_args(script_args)
	return parsed_script_args


def mk_file_path(orig_fname, suffix):
	parts = orig_fname.split(".")
	if len(parts) < 2:
		return orig_fname + "-" + suffix + ".stl"
	return parts[0] + "-" + suffix + "." + parts[1]

class Piece:
	def __init__(self):
		self.start_height = float(args.start_height)
		if args.end_height:
			self.end_height = float(args.end_height)
		else:
			self.end_height = float("inf")

	def get_fname(self):
		return self.get_piece_name().lower() + ".stl"

	def get_edge_w(self, x):
		return edge_w

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

	def make_prismic_part(self, start_x, part_h, vn_arg, cur_poly_f, part_name, all_in_one = False):
		if not all_in_one:
			mu.reset_scene()
		n_poly = None
		x = start_x
		dx = part_h / float(vn_arg)
		verts = []
		faces = []

		for v_ind in range(0, vn_arg + 1):
			cur_poly = cur_poly_f(x)
			verts += cur_poly
			x += dx
			if n_poly == None:
				n_poly = len(cur_poly)
			if v_ind > 0:
				base_low = (v_ind - 1) * n_poly
				for i in range(0, n_poly):
					next_i = base_low + (i + 1) % n_poly
					faces.append([ base_low + i, next_i, next_i + n_poly,
								base_low + i + n_poly ])
		faces.append(range(0, n_poly))
		faces.append(range((vn_arg) * n_poly, (vn_arg + 1) * n_poly))
		mu.create_mesh_from_data(self.get_piece_name() + ' ' + part_name, o_v, verts, faces, True)
		if not all_in_one:
			self.save_part(part_name)

	def make_conic_part(self, start_x, part_h, vn_arg, cur_r_f, part_name, all_in_one = False):
		if not all_in_one:
			mu.reset_scene()
		self.make_conic_part_low(start_x, part_h, vn_arg,
															 cur_r_f, part_name)
		if not all_in_one:
			self.save_part(part_name)

	def out_of_bounds(self, x):
		# print("x="+str(x)+", start_height=" + str(self.start_height))
		return x < self.start_height or x > self.end_height

	def make_conic_part_low(self, start_x, part_h, vn_arg, cur_r_f, part_name):
		x = start_x
		dx = part_h / float(vn_arg)
		top_x = start_x + part_h
		verts = []
		inner_verts = []
		faces = []
		inner_faces = []
		outer_faces = []
		top_inner_face = None
		top_outer_face = None
		bottom_inner_face = None
		bottom_outer_face = None
		cur_r = None
		first_r = False
		skipped_rings = 0
		print("part_name = " + part_name + " start_x="+str(x) +
			" dx=" + str(dx) + " start_height=" + str(self.start_height))

		for v_ind in range(0, vn_arg + 1):
			# out_of_bounds normally should always return false,
			# but we might have a case when the print completed partially
			# and we want to print the remaining part
			if self.out_of_bounds(x):
				x += dx
				skipped_rings += 1
				continue
			if cur_r == None:
				first_r = True
			cur_r = cur_r_f(x)
			if first_r:
				print("start_r = " + str(cur_r))
				first_r = False
			if cur_r <= 0:
				break
			# print("x=" + str(x) + ", cur_r = " + str(cur_r))
			cur_verts = mu.get_circle_verts((0,0,x), cur_r, n)
			cur_n_verts = len(verts)
			verts += cur_verts
			inner_r = cur_r - self.get_edge_w(x)
			if inner_r < 0:
				inner_r = 0.0
			cur_inner_verts = mu.get_circle_verts((0,0,x), inner_r, n)
			inner_verts += cur_inner_verts
			x += dx
			real_v_ind = v_ind - skipped_rings

			if real_v_ind > 0:
				base_low = (real_v_ind - 1) * n
				for i in range(0, n):
					next_i = base_low + (i + 1) % n
					#print(base_low)
					face = [ base_low + i, next_i, next_i + n,
								base_low + i + n ]
					inner_face = face.copy()
					inner_face.reverse()
					faces.append(face)
					inner_faces.append(inner_face)

		top_outer_face = list(range(0,n))
		for i,f in enumerate(inner_faces):
			for j,v in enumerate(f):
				inner_faces[i][j] += len(verts)
		top_inner_face = list(range(len(verts), len(verts) + n))
		bottom_outer_face = list(range(len(verts) - n, len(verts)))
		bottom_inner_face = list(range(2 * len(verts) - n, 2 * len(verts)))
		verts += inner_verts
		faces += inner_faces


		if bottom_outer_face:
			faces += mu.connect_circles(bottom_outer_face, bottom_inner_face, 0)
		if top_outer_face:
			faces += mu.connect_circles(top_outer_face, top_inner_face, 0)

		print("final x="+str(x) + " final r = " + str(cur_r))
		print("n_verts="+str(len(verts)) + " skipped_rings = " + str(skipped_rings))
		#print("faces=" + str(faces))
		#print("inner_faces=" + str(inner_faces))

		print("n_faces="+str(len(faces)))
		mu.create_mesh_from_data(self.get_piece_name() + ' ' + part_name, o_v, verts, faces, True)

# three-part part piece
class Piece3(Piece):
	def __init__(self):
		super().__init__()
		self.base_h = h * self.base_q()
		self.part = args.part
		self.middle_h = h * self.middle_q()
		self.top_h = h * self.total_q() - self.base_h - self.middle_h
		self.base_top_r = r * self.top_base_rq()
		self.bend_h = self.base_h + self.middle_h * self.bend_q()
		self.bend_r = self.base_top_r * self.bend_rq()
		self.vn_base = int(vn * self.base_q())
		self.vn_middle = int(vn * self.middle_q())
		self.vn_top = int(vn * self.top_q())
		self.middle_top_r = self.middle_cur_r(self.middle_h + self.base_h)

	def top_base_rq(self):
		return 0.5
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
			for p in ["base", "middle", "top"]:
				if self.part != "all" and self.part != p:
					continue
				getattr(self, "make_" + p)(all_in_one)
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
	def top_x_fix(self):
		return self.top_h * 0.0

	# top_shift_h + top_r = top_h
	# (top_shift_h + top_r) = 2 * top_r * top_ball_q
	# top_h = 2 * top_r * top_ball_q
	# top_r = top_h /( 2 * top_ball_q)
	def top_cur_r(self, x):
		d_fix = self.top_x_fix()
		dx = x - self.top_shift_h - self.middle_h - self.base_h + d_fix
		r_fix = self.top_r - d_fix
		y_sq = r_fix * r_fix - dx * dx
		#print("y_sq = " + str(y_sq) + " x = " + str(x))
		if y_sq > 0.0:
			return math.sqrt(y_sq)
		else:
			return 0

# model knight:
# base_r = 9.23
# base_h = 12.42
# support_r = 5.62
# support_h = 2.2
# bottom_neck_w = 13.1
# bottom_neck_l = 15.1
# bend_neck_w = 10.4
# bend_neck_l = 16.1
class Knight(Piece3):
	def __init__(self):
		super().__init__()
		self.base_poly_d = self.base_top_r / math.sqrt(2)
		self.poly_n = n
		self.support_h = self.middle_h * self.support_q()
		self.support_rq = 0.7
		self.support_r = self.support_rq * self.base_top_r
		self.head_q = 1.0
		self.bend_neck_w = self.support_r * self.base_support_q() * (1 -
				self.middle_low_yq())
		self.top_neck_w = self.bend_neck_w * self.top_neck_q()
		self.top_poly_l = self.base_poly_d * self.top_dq()
		self.top_poly_w = self.top_neck_w
		self.middle_bend_q = 0.18
		self.middle_bend_h = self.middle_h * self.middle_bend_q
		self.bend_h2 = self.middle_h - self.middle_bend_h
		self.bend_neck_l = self.support_r * self.middle_upper_q()
		print("bend_neck_l=" + str(self.bend_neck_l))
		self.middle_ell_y = self.middle_bend_h * ( 1 + self.middle_circle_shift_q())
		self.ear_h = self.top_h * self.ear_q()
		self.ear_l = self.top_poly_l * self.ear_ql()
		self.ear_w = self.top_poly_w * self.ear_qw()
		self.eye_h = self.top_h * self.eye_q()
		self.eye_h_pos = self.top_h - 2 * self.eye_h - self.ear_h
		self.eye_l_pos = 2 * self.top_poly_l - self.ear_l * 4
		self.eye_w = self.eye_h
		self.eye_l = self.eye_h * 2

	def get_piece_name(self):
		return "Knight"
	def eye_q(self):
		return 0.1
	def eye_ql(self):
		return 0.1
	def eye_qw(self):
		return 0.1
	def top_neck_q(self):
		return 0.6
	def top_dwq(self):
		return 0.5
	def top_dq(self):
		return 1.2
	def support_q(self):
		return 0.1
	def top_base_rq(self):
		return 0.7
	def total_q(self):
		return 1.1
	def base_q(self):
		return 0.2
	def middle_q(self):
		return 0.65
	def middle_x_back_q(self):
		return 2.0
	def middle_circle_shift_q(self):
		return 0.2
	def middle_low_yq(self):
		return 0.5
	def base_support_q(self):
		return 1.1
	def middle_upper_q(self):
		return 1.1
	def bend_neck_q(self):
		return 2.0

	def middle_cur_poly(self, z):
		if z < self.base_h + self.support_h:
			return mu.get_circle_verts((0, 0, z), self.support_r, self.poly_n)
		dz = z - self.base_h - self.support_h
		#dx = self.base_poly_d * dz/self.middle_h
		if dz < self.middle_bend_h:
			dx = self.bend_neck_l * math.sqrt(1 -
				((self.middle_bend_h - dz)/(self.middle_ell_y)) ** 2)
			dx_back = -dx
			dy = self.support_r * self.base_support_q() * (1 - dz / self.middle_bend_h *
				self.middle_low_yq())
		else:
			dx_back = -self.bend_neck_l
			dx =  (dz - self.middle_bend_h) / (self.middle_h - self.middle_bend_h) * (1 -
				self.bend_neck_q()) * self.bend_neck_l + self.bend_neck_l
			#print("dx=" + str(dx))
			dy = (self.top_neck_w - self.bend_neck_w)* ((dz - self.middle_bend_h)/
				(self.middle_h - self.middle_bend_h)) + self.bend_neck_w
		return mu.extend_poly([
														[dx_back, -dy, z],
														[dx_back, dy, z],
														[dx, dy, z] ,
														[dx, -dy, z]
													], self.poly_n)

	def top_shift_q(self):
		return 0.2
	def ear_q(self):
		return 0.4
	def ear_qw(self):
		return 0.2
	def ear_ql(self):
		return 0.2

	def ear_poly_left(self, z):
		return self.ear_poly_low(z, -self.top_poly_l + 2 * self.ear_l,
														 self.top_poly_w - 2 * self.ear_w)

	def ear_poly_right(self, z):
		return self.ear_poly_low(z, -self.top_poly_l + 2 * self.ear_l,
														-self.top_poly_w + 2 * self.ear_w)

	def ear_poly_low(self, z, x, y):
		dz = self.top_h + self.base_h + self.middle_h - z
		k = dz / self.ear_h
		dx = self.ear_l * k
		dy = self.ear_w * k
		poly = [
			[x + dx, y + dy, z], [x - dx, y + dy, z],
			[x - dx, y - dy, z], [x + dy, y - dy, z]
		]
		return mu.extend_poly(poly, self.poly_n)

	def add_eyes_to_poly_low(self, poly, z, shift_x, cur_eye_w, cur_eye_l, cur_eye_dl):
		eye_1_poly = [
			[self.top_poly_l + shift_x - self.eye_l_pos, -self.top_poly_w, z],
			[self.top_poly_l + shift_x - self.eye_l_pos - cur_eye_dl, -self.top_poly_w + cur_eye_w, z],
			[self.top_poly_l + shift_x - self.eye_l_pos - self.eye_l + cur_eye_dl, -self.top_poly_w + cur_eye_w, z],
			[self.top_poly_l + shift_x - self.eye_l_pos - self.eye_l, -self.top_poly_w, z]
		]
		eye_2_poly = [
			[self.top_poly_l + shift_x - self.eye_l_pos - self.eye_l, self.top_poly_w, z],
			[self.top_poly_l + shift_x - self.eye_l_pos - self.eye_l + cur_eye_dl,
					self.top_poly_w - cur_eye_w, z],
			[self.top_poly_l + shift_x - self.eye_l_pos - cur_eye_dl , self.top_poly_w - cur_eye_w, z],
			[self.top_poly_l + shift_x - self.eye_l_pos , self.top_poly_w, z]
		]
		return poly[0:2] + eye_1_poly + poly[2:] + eye_2_poly + poly[0:1]

	def add_eyes_to_poly(self, poly, z, shift_x):
		eye_r = self.eye_h/2
		rel_z = z - self.middle_h - self.base_h - self.eye_h_pos
		q = (1 - math.sqrt(1 - ((rel_z - eye_r)/eye_r)**2))
		cur_eye_w = self.eye_w * q
		cur_eye_l = self.eye_l * (1 - q)
		cur_eye_dl = (self.eye_l - cur_eye_l)/2
		return self.add_eyes_to_poly_low(poly, z, shift_x, cur_eye_w, cur_eye_l, cur_eye_dl)

	def top_cur_poly(self, z):
		dx = 0
		#shift_x = self.top_poly_l * self.top_shift_q()
		shift_x = 0
		dz = z - self.middle_h - self.base_h
		poly = [[self.top_poly_l + shift_x, self.top_poly_w, z],
					[self.top_poly_l + shift_x, -self.top_poly_w, z],
					[-self.top_poly_l + dx + shift_x, -self.top_poly_w , z] ,
					[-self.top_poly_l + dx + shift_x, self.top_poly_w , z]
				]
		print("dz=" + str(dz) + "eye_h_pos=" + str(self.eye_h_pos) + " top_h")
		if dz > self.eye_h_pos and dz < self.eye_h_pos + self.eye_h:
			poly = self.add_eyes_to_poly(poly, z, shift_x)
			print("eye_poly: " + str(poly))
		else:
			# a hack to align the polygons at the joint with eyes
			poly = self.add_eyes_to_poly_low(poly, z, shift_x, 0, self.eye_l, self.eye_l)
		poly = mu.extend_poly(poly, self.poly_n)
		#print("poly has length " + str(len(poly)))
		#print("poly = " + str(poly))
		return poly

	def make_top(self, all_in_one):
		if not all_in_one:
			mu.reset_scene()
		self.make_prismic_part(self.middle_h + self.base_h, self.top_h - self.ear_h, self.vn_top,
												self.top_cur_poly, "top", True)
		for f in [self.ear_poly_left, self.ear_poly_right]:
			self.make_prismic_part(self.middle_h + self.base_h + self.top_h - self.ear_h,
													self.ear_h, self.vn_top, f, "top", True)
		if not all_in_one:
			self.save_part("top")


	def make_middle(self, all_in_one):
		self.make_prismic_part(self.base_h, self.middle_h, self.vn_middle,
												self.middle_cur_poly, "middle", all_in_one)


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
	def get_edge_w(self, x):
		save_x = x
		x -= self.base_h + self.middle_h
		if x < 0:
			return edge_w
		#print("funnel_h = " + str(self.funnel_h) + ", x = " + str(x))
		cur_r = self.top_cur_r(save_x)
		#print("funnel_h = " + str(self.funnel_h) + ", x = " + str(x) + " cur_r = " + str(cur_r))
		if self.top_h - x < edge_w:
			#print(" returning cur_r = " + str(cur_r))
			return cur_r
		res =  cur_r - self.crown_r * (self.top_h - x)/self.top_h
		#print("returning " + str(res))
		if res < edge_w:
			return edge_w
		return res

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
		r = self.crown_r * (self.top_h - x)/(self.top_h - self.funnel_h)
		return r

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
		self.dent_h = self.top_h / self.vn_top # indention to make assembly easier
		self.middle_h += self.dent_h
		print("Rook: pit_h = " + str(self.pit_h) + " pit_r = " + str(self.pit_r) +
			" top_r = " + str(self.top_r) + " dent_h = " + str(self.dent_h)
		)
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
	def get_edge_w(self, x):
		x_rel_top = x - self.middle_h - self.base_h + self.dent_h
		# we are below the top part
		#print("x_rel_top=" + str(x_rel_top))
		if x_rel_top <= self.dent_h or x_rel_top > self.pit_h:
			return edge_w
		return self.top_cur_r(x)

def bail(msg):
	print("Error: " + msg)
	sys.exit(1)

def get_piece_obj(name):
	name = name[0].upper() + name[1:].lower()
	try:
		obj = globals()[name]
	except:
		bail("Bad piece name " + name)
	return obj()

def make_all():
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

args = get_args()

origin = (0,0,0)
o_v = Vector(origin)

h = float(args.height)
r = float(args.radius)
th = float(args.thickness)
n = int(args.num_parts)
vn = int(args.num_vparts)
edge_w = float(args.edge_w)
piece = args.piece

if piece:
	piece_obj = get_piece_obj(piece)
	piece_obj.make()
else:
	make_all()
