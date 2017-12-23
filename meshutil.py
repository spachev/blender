import math
import bpy
from math import radians
from mathutils import Matrix,Vector
from itertools import product

def mul_v(v,k):
	return [[sub_el * k for sub_el in el] for el in v]

def get_brick_faces():
	return [[0,1,3,2], [6,7,5,4],
		[5,7,3,1], [4,0,2,6],
		[0,4,5,1],[6,2,3,7]
	]

def get_brick_verts(o,l,w,h):
	verts = []
	plus_minus = (-1.0,1.0)

	for x,y,z in product(plus_minus, repeat=3):
		verts.append((x*l/2.0, y*w/2.0, z*h/2.0))

	return verts

def get_hollow_brick_faces():
	return [[5,7,3,1], [4,0,2,6],
		[0,4,5,1],[6,2,3,7]
	]

# dl - delta length for the top part of the trapezoidal prism
def get_tp_prism_verts(o,l,w,h,dl,rot_x,rot_y,rot_z):
	verts = []
	plus_minus = (-1.0,1.0)

	for x,y,z in product(plus_minus, repeat=3):
		if x < 0.0 and z > 0.0:
			verts.append(Vector((x*l/2.0 , y*w/2.0, z*h/2.0)))
		else:
			verts.append(Vector((x*(l/2.0+dl), y*w/2.0, z*h/2.0)))

	mat_rot_x = Matrix.Rotation(rot_x, 3, 'X')
	mat_rot_y = Matrix.Rotation(rot_y, 3, 'Y')
	mat_rot_z = Matrix.Rotation(rot_z, 3, 'Z')
	mat_rot = mat_rot_x * mat_rot_y * mat_rot_z
	print(mat_rot)
	return [v * mat_rot  for v in verts]

def get_circle_verts(o,r,n):
	theta = 0
	verts = []
	pi2 = 2.0 * math.pi
	d_theta = pi2 / n

	while theta < pi2:
		x = o[0] + r * math.cos(theta)
		y = o[1] + r * math.sin(theta)
		verts.append([x,y,o[2]])
		theta += d_theta

	return verts

def extend_poly(poly, target_n):
	n_poly = len(poly)
	if n_poly > target_n:
		raise Exception("Polygon has too many verts")
	verts_per_side = target_n // n_poly
	if target_n % n_poly:
		verts_per_side += 1
	new_poly = []
	for i in range(0, n_poly):
		v1 = poly[i]
		v2 = poly[(i + 1) % n_poly]
		dv = [0, 0, 0]
		for j in range(0,3):
			dv[j] = (v2[j] - v1[j]) / verts_per_side
		for k in range(0, verts_per_side):
			v = [0, 0, 0]
			for j in range(0,3):
				v[j] = v1[j] + k * dv[j]
			new_poly.append(v)
			if len(new_poly) == target_n:
				return new_poly
	raise Exception("Bug: should never get here: len(new_poly) = " + str(len(new_poly)))


def reset_scene():
	scn = bpy.context.scene
	obs = scn.objects
	for ob in obs:
		obs.unlink(ob)

def create_mesh_from_data(name, origin, verts, faces, remove_cube=True):
	me = bpy.data.meshes.new(name+'Mesh')
	ob = bpy.data.objects.new(name, me)
	ob.location = origin
	ob.show_name = True

	# Link object to scene and make active
	scn = bpy.context.scene
	obs = scn.objects

	if remove_cube:

		cube_obs = [cube_ob for cube_ob in obs if cube_ob.name == 'Cube']
		if len(cube_obs):
			obs.unlink(cube_obs[0])

	scn.objects.link(ob)
	scn.objects.active = ob
	ob.select = True

	# Create mesh from given verts, faces.
	me.from_pydata(verts, [], faces)
	# Update mesh with new data
	me.update()
	return me

def remove_cube():
	scn = bpy.context.scene
	obs = scn.objects
	cube_ob = [cube_ob for cube_ob in obs if cube_ob.name == 'Cube'][0]
	if cube_ob:
		obs.unlink(cube_ob)

def add_text_new(txt, origin, r, n, font_size):
	bpy.ops.object.text_add(
		location=origin,
		rotation=(math.pi/4,0,0))
	ob = bpy.context.object
	ob.name = 'NameText'
	tcu = ob.data
	tcu.name = 'NameTextData'
	pl = tcu.splines.new('BEZIER')
	circle = get_circle_verts(origin, r, n)
	pl.bezier_points.add(n)

	for i,p in enumerate(circle):
		pl.bezier_points[i].co = (p[0], p[1], p[2])

	print(dir(tcu))
	# TextCurve attributes
	tcu.body = txt
	tcu.font = bpy.data.fonts[0]
	tcu.offset_x = -9
	tcu.offset_y = -0.25
	tcu.shear = 0.5
	tcu.size = font_size
	tcu.space_character = 2
	tcu.space_word = 4
	tcu.use_path = True
	print(tcu.use_path)

	print(tcu.splines)
	# Inherited Curve attributes
	tcu.extrude = font_size/10.0
	tcu.fill_mode="FRONT"
	tcu.use_fill_deform = True
	tcu.fill_mode="FRONT"

def create_bevel_object(coords):
	# Create Bevel curve and object
	cu = bpy.data.curves.new('BevelCurve', 'CURVE')
	ob = bpy.data.objects.new('BevelObject', cu)
	bpy.context.scene.objects.link(ob)

	# Set some attributes
	cu.dimensions = '3D'
	cu.resolution_u = 6
	cu.twist_mode = 'MINIMUM'
	ob.show_name = True


	# Create spline and set control points
	spline = cu.splines.new('NURBS')
	nPointsU = len(coords)
	spline.points.add(nPointsU)
	for n in range(nPointsU):
		spline.points[n].co = coords[n] + [1.0,]

	# Set spline attributes. Points probably need to exist here.
	spline.use_cyclic_u = True
	spline.resolution_u = 1
	spline.order_u = 3

	return ob

def get_curve(coords):
	curveData = bpy.data.curves.new('myCurve', type='CURVE')
	curveData.dimensions = '3D'
	curveData.resolution_u = 2

	# map coords to spline
	polyline = curveData.splines.new('POLY')
	polyline.points.add(len(coords)-1)
	for i, coord in enumerate(coords):
		x,y,z = coord
		polyline.points[i].co = (x, y, z, 1)

	curveOB = bpy.data.objects.new('myCurve', curveData)
	# attach to scene and validate context
	scn = bpy.context.scene
	scn.objects.link(curveOB)
	scn.objects.active = curveOB
	curveOB.select = True
	return curveOB

def add_text(txt, origin, r, n, font_size):
	scene = bpy.context.scene
	mat_rot_x = Matrix.Rotation(radians(90.0), 4, 'X')
	mat_rot_y = Matrix.Rotation(radians(0.0), 4, 'Y')
	fc = bpy.data.curves.new(type="FONT", name="fontCurve") # todo - figure out name
	#print(dir(fc))
	circle = get_circle_verts((0,0,0), r, n)
	fo = bpy.data.objects.new("myFontOb", fc)
	#print("Adding " + txt)
	#print(fo.data)
	fo.data.body = txt
	fo.data.size = font_size
	fo.data.extrude = font_size/10.0
	fo.data.bevel_depth = 0.1
	scene.objects.link(fo)
	scene.objects.active = fo
	fo.select = True
	fo.location = origin
	fo.rotation_euler=(radians(90),radians(0),radians(0))
	m = fo.modifiers.new('curve', 'CURVE')
	m.object = get_curve(circle)
	#print(fo.modifiers)
	#fo.data.follow_curve = get_curve(circle)
	#fo.matrix_world = mat_rot_x * mat_rot_y
	scene.update()
	#print(len(scene.objects))

