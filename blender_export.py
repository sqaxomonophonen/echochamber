import bpy
from bpy_extras import mesh_utils
from mathutils import Vector, Matrix
import sys
import os
import math
from struct import pack

src = bpy.data.filepath
dst = sys.argv[sys.argv.index('--') + 1]

print("exporting %s => %s ..." % (src, dst))


class ISet:
	def __init__(self):
		self.vimap = {}
		self.lst = []

	def _add(self, v):
		if v not in self.vimap:
			self.vimap[v] = len(self.lst)
			self.lst.append(v)
		return self.vimap[v]

	def write_block(self, b):
		b.u32(len(self.lst))
		for data in self.lst:
			self._write_data(b, data)

class Biquads(ISet):
	def get_block_id(self):
		return b"BQD0"

	# http://www.earlevel.com/main/2003/02/28/biquads/
	# http://www.earlevel.com/main/2013/10/13/biquad-calculator-v2/
	def add(self, a0, a1, a2, b1, b2):
		return self._add((a0, a1, a2, b1, b2))

	def add_std_highpass(self):
		# hz=48000hz
		# fc=1375hz
		# q=0.7071
		a0 = 0.8804686573973374
		a1 = -1.7609373147946747
		a2 = 0.8804686573973374
		b1 = -1.7465984140629678
		b2 = 0.7752762155263814
		# flipping a0/a1/a2 to cause phase inversion. THIS IS HIGHLY
		# SCIENTIFIC I SAW IT ON FRINGE
		return self.add(-a0, -a1, -a2, b1, b2)

	def _write_data(self, b, data):
		a0,a1,a2,b1,b2 = data
		b.f32(a0)
		b.f32(a1)
		b.f32(a2)
		b.f32(b1)
		b.f32(b2)


class EmissionGroups(ISet):
	def get_block_id(self):
		return b"EMG0"

	def add(self, name = None):
		return self._add(name or "emit000")

	def _write_data(self, b, data):
		b.padded_string(data, 64)

class Microphones:
	def get_block_id(self):
		return b"MIC0"

	def __init__(self):
		self.lst = []
		self.name_set = set()

	def add(self, position, name = None):
		if name is None:
			name = "mic%.3d" % len(self.lst)
		if name in self.name_set:
			raise RuntimeError("camera name conflict; '%s' already used" % name)
		self.name_set.add(name)
		self.lst.append((name, position))

	def write_block(self, b):
		b.u32(len(self.lst))
		for (name, position) in self.lst:
			b.padded_string(name, 64)
			b.vec3(position)

class Materials(ISet):
	def get_block_id(self):
		return b"MAT0"

	def add(self, emission_group_id, biquad_id, hardness):
		return self._add((emission_group_id, biquad_id, hardness))

	def _write_data(self, b, data):
		emission_group_id, biquad_id, hardness = data
		b.i32(emission_group_id)
		b.u32(biquad_id)
		b.f32(hardness)

class Polys:
	def get_block_id(self):
		return b"PLY0"

	def __init__(self):
		self.lst = []

	def add(self, material_id, coords):
		self.lst.append((material_id, coords))

	def write_block(self, b):
		for material_id, coords in self.lst:
			b.u32(material_id)
			b.u32(len(coords))
			for c in coords:
				b.vec3(c)

biquads = Biquads()
emission_groups = EmissionGroups()
microphones = Microphones()
materials = Materials()
polys = Polys()


for bo in bpy.data.objects:
	if bo.type == "MESH":
		mesh = bo.to_mesh(bpy.context.scene, True, "PREVIEW")
		for polygon in mesh.polygons:
			coords = []
			for i in range(polygon.loop_start, polygon.loop_start + polygon.loop_total):
				coord = bo.matrix_world * mesh.vertices[mesh.loops[i].vertex_index].co
				coords.append(coord)

			# flip polygon if normal is wrong
			e0 = coords[1] - coords[0]
			e1 = coords[1] - coords[2]
			if polygon.normal.dot(e1.cross(e0).normalized()) < 0:
				coords.reverse()

			material = mesh.materials[polygon.material_index]

			hardness = material.specular_hardness/512.0

			if material.emit > 0:
				emission_group_id = emission_groups.add()
				biquad_id = biquads.add(material.emit, 0.0, 0.0, 0.0, 0.0)
			else:
				emission_group_id = -1
				# TODO use material.diffuse_color to create a
				# "default" biquad in the future?
				# TODO also allow direct definition via custom
				# properties
				biquad_id = biquads.add_std_highpass()

			material_id = materials.add(emission_group_id, biquad_id, hardness)

			polys.add(material_id, coords)

	elif bo.type == "CAMERA":
		loc, rot, scale = bo.matrix_world.decompose()
		microphones.add(loc)


class ECSBlock:
	def __init__(self):
		self.buf = b""
		self.id = None

	def set_id(self, id):
		self.id = id

	def u32(self, i):
		self.buf += pack("I", i)

	def i32(self, i):
		self.buf += pack("i", i)

	def f32(self, i):
		self.buf += pack("f", i)

	def vec3(self, v):
		for f in v:
			self.f32(f)

	def padded_string(self, s, padding):
		if len(s) >= padding:
			raise RuntimeError("string %s exceeds padding length" % s)
		self.buf += bytes(s, "ascii")
		self.buf += b"\x00" * (padding - len(s))

	def pack(self):
		if not isinstance(self.id, bytes) or len(self.id) != 4:
			raise RuntimeError("invalid block id %s" % self.id)
		o = self.id
		o += pack("I", len(self.buf))
		o += self.buf
		o += pack("I", 0xdeadbeef)
		return o



class ECS:
	def __init__(self, filename):
		self.filename = filename

	def write(self, objs):
		# header
		buf = b"ECSn"
		buf += pack("I", 0xec511235) # endianess detect
		buf += pack("I", 1) # version

		# blocks
		for obj in objs:
			block = ECSBlock()
			block.set_id(obj.get_block_id())
			obj.write_block(block)
			buf += block.pack()

		with open(self.filename, "wb") as f:
			f.write(buf)

ECS(dst).write([biquads, emission_groups, microphones, materials, polys])
