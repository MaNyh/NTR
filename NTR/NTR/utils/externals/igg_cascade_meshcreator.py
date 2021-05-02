import sys
import os
import pickle

print("starting")
print(__file__)

def read_pickle_args(path):
    filepath = os.path.join(path,"args.pkl")
    with open(filepath,"rb") as Fobj:
        dict = pickle.load(Fobj)
    return dict


script_path = os.path.dirname(os.path.abspath(__file__))
args = read_pickle_args(script_path)

pointcloudfile = args["pointcloudfile"]
add_path = args["add_path"]
sys.path.append(add_path)

from utils.externals.tecplot_functions import openTecplotFile


yPerLowHGridBlockPitchStart = args["yPerLowHGridBlockPitchStart"]
yPerHighHGridBlockPitchStart = args["yPerHighHGridBlockPitchStart"]
vk_BlockStartFromChord = args["vk_BlockStartFromChord"]
hk_BlockStartFromChord = args["hk_BlockStartFromChord"]
factor = args["factor"]

delta_i = args["delta_i"]
delta_i /= factor

cellwidthcoeff = args["cellwidthcoeff"]

first_cell_width = delta_i*args["first_cell_width"]
first_cell_width *= factor

exp_ratio = args["exp_ratio"]

layers = args["layers"]
layers = int(layers*factor)

extrudeLength = args["extrudeLength"]
extrudeNodes = args["extrudeNodes"]
extrudeNodes *= int(factor*extrudeNodes)

print("factor: " +str(factor))
print("delta_i: " +str(delta_i))
print("first_cell_width: " +str(first_cell_width))
print("exp_ratio: " +str(exp_ratio))
print("layers: " +str(layers))
print("extrudeLength: " +str(extrudeLength))
print("extrudeNodes: " +str(extrudeNodes))


data = openTecplotFile(pointcloudfile)
x_ss, y_ss = data[0][0], data[0][1]
x_ps, y_ps = data[1][0], data[1][1]
x_lower, y_lower = data[2][0], data[2][1]
x_upper, y_upper = data[3][0], data[3][1]
x_sc, y_sc = data[4][0], data[4][1]

pitch = y_upper[0] - y_lower[0]


cspline_peri_lower = new_cspline("cspline_peri_lower")
cspline_peri_upper = new_cspline("cspline_peri_upper")
cspline_ss = new_cspline("cspline_ss")
cspline_ps = new_cspline("cspline_ps")
cspline_inlet = new_cspline("cspline_inlet")
cspline_outlet = new_cspline("cspline_outlet")

for i in range(len(x_lower)):
    cspline_peri_lower.insert_point(i + 1, Point(x_lower[i], y_lower[i], 0))
    cspline_peri_upper.insert_point(i + 1, Point(x_upper[i], y_upper[i], 0))

for i in range(len(x_ss)):
    cspline_ss.insert_point(i + 1, Point(x_ss[i], y_ss[i], 0))

for i in range(len(x_ps)):
    cspline_ps.insert_point(i + 1, Point(x_ps[i], y_ps[i], 0))

cspline_ps.insert_point(len(x_ps) + 1, Point(x_ss[0], y_ss[0], 0))
cspline_inlet.insert_point(1, Point(x_lower[0], y_lower[0], 0))
cspline_inlet.insert_point(2, Point(x_upper[0], y_upper[0], 0))
cspline_outlet.insert_point(1, Point(x_lower[-1], y_lower[-1], 0))
cspline_outlet.insert_point(2, Point(x_upper[-1], y_upper[-1], 0))

# =============================================================================
# Blockerstellung
# =============================================================================


# p1=CurvePointNorm(Curve("cspline_peri_upper"),Curve("cspline_peri_upper").calc_normalize(Curve("cspline_peri_upper").project_point(CurvePointNorm(Curve("cspline_ps"),0.0))))
p1 = CurvePointNorm(Curve("cspline_peri_upper"),
                    Curve("cspline_peri_upper").calc_normalize(
                        Curve("cspline_peri_upper").project_point(Point(CurvePointNorm(Curve("cspline_ps"), 0.0).x,
                                                                        CurvePointNorm(Curve("cspline_ps"),
                                                                                       0.0).y + 0.4 * pitch,
                                                                        0))))



"""
p2=CurvePointNorm(Curve("cspline_peri_lower"),
		  	Curve("cspline_peri_lower").calc_normalize(Curve("cspline_peri_lower").project_point(CurvePointNorm(Curve("cspline_ss"),
			0.0))))
"""



p2 = CurvePointNorm(Curve("cspline_peri_lower"),
                    Curve("cspline_peri_lower").calc_normalize(
                        Curve("cspline_peri_lower").project_point(Point(CurvePointNorm(Curve("cspline_ps"), 1.0).x,
                                                                        CurvePointNorm(Curve("cspline_ps"),
                                                                                       0.0).y + 0.15 * pitch,
                                                                        1.0))) + 0.05)

p_vk = CurvePointNorm(Curve("cspline_ps"), 0.0)
p_hk = CurvePointNorm(Curve("cspline_ps"), 1.0)

new_block_face(CurvePointNorm(Curve("cspline_peri_lower"), 0.0),
               CurvePointNorm(Curve("cspline_peri_lower"), Curve("cspline_peri_lower").calc_normalize(
                   Curve("cspline_peri_lower").project_point(Point(p1.x, p1.y - pitch, 0)))),
               Point(p1.x, p1.y + yPerLowHGridBlockPitchStart * pitch - pitch, 0),
               CurvePointNorm(Curve("cspline_inlet"), yPerLowHGridBlockPitchStart)
               )

new_block_face(CurvePointNorm(Curve("cspline_inlet"), yPerLowHGridBlockPitchStart),
               Point(p1.x, p1.y + yPerLowHGridBlockPitchStart * pitch - pitch, 0),
               Point(p1.x, p1.y + yPerHighHGridBlockPitchStart * pitch - pitch, 0),
               CurvePointNorm(Curve("cspline_inlet"), yPerHighHGridBlockPitchStart)
               )

new_block_face(CurvePointNorm(Curve("cspline_inlet"), yPerHighHGridBlockPitchStart),
               Point(p1.x, p1.y + yPerHighHGridBlockPitchStart * pitch - pitch, 0),
               p1,
               CurvePointNorm(Curve("cspline_peri_upper"), 0.0)
               )

new_block_face(CurvePointNorm(Curve("cspline_peri_lower"), Curve("cspline_peri_lower").calc_normalize(
    Curve("cspline_peri_lower").project_point(Point(p1.x, p1.y - pitch, 0)))),
               CurvePointNorm(Curve("cspline_peri_lower"), Curve("cspline_peri_lower").calc_normalize(
                   Curve("cspline_peri_lower").project_point(Point(p2.x, p2.y, 0)))),
               Point(p2.x, p2.y + yPerLowHGridBlockPitchStart * pitch, 0),
               Point(p1.x, p1.y + yPerLowHGridBlockPitchStart * pitch - pitch, 0)
               )

new_block_face(CurvePointNorm(Curve("cspline_peri_lower"), Curve("cspline_peri_lower").calc_normalize(
    Curve("cspline_peri_lower").project_point(Point(p2.x, p2.y, 0)))),
               CurvePointNorm(Curve("cspline_peri_lower"), 1.0),
               CurvePointNorm(Curve("cspline_outlet"), yPerLowHGridBlockPitchStart),
               Point(p2.x, p2.y + yPerLowHGridBlockPitchStart * pitch, 0)
               )

new_block_face(Point(p2.x, p2.y + yPerLowHGridBlockPitchStart * pitch, 0),
               CurvePointNorm(Curve("cspline_outlet"), yPerLowHGridBlockPitchStart),
               CurvePointNorm(Curve("cspline_outlet"), yPerHighHGridBlockPitchStart),
               Point(p2.x, p2.y + yPerHighHGridBlockPitchStart * pitch, 0)
               )

new_block_face(Point(p2.x, p2.y + yPerHighHGridBlockPitchStart * pitch, 0),
               CurvePointNorm(Curve("cspline_outlet"), yPerHighHGridBlockPitchStart),
               CurvePointNorm(Curve("cspline_peri_upper"), 1.0),
               CurvePointNorm(Curve("cspline_peri_upper"), Curve("cspline_peri_upper").calc_normalize(
                   Curve("cspline_peri_upper").project_point(Point(p2.x, p2.y + pitch, 0))))
               )

new_block_face(Point(p1.x, p1.y + yPerHighHGridBlockPitchStart * pitch - pitch, 0),
               Point(p2.x, p2.y + yPerHighHGridBlockPitchStart * pitch, 0),
               CurvePointNorm(Curve("cspline_peri_upper"), Curve("cspline_peri_upper").calc_normalize(
                   Curve("cspline_peri_upper").project_point(Point(p2.x, p2.y + pitch, 0)))),
               p1
               )

# O-Netz
# Vk

new_block_face(Point(p1.x, p1.y + yPerLowHGridBlockPitchStart * pitch - pitch, 0),
               CurvePointNorm(Curve("cspline_ps"), vk_BlockStartFromChord),
               CurvePointNorm(Curve("cspline_ss"), hk_BlockStartFromChord),
               Point(p1.x, p1.y + yPerHighHGridBlockPitchStart * pitch - pitch, 0)
               )

# Hk

new_block_face(Point(p2.x, p2.y + yPerLowHGridBlockPitchStart * pitch, 0),
               Point(p2.x, p2.y + yPerHighHGridBlockPitchStart * pitch, 0),
               CurvePointNorm(Curve("cspline_ss"), vk_BlockStartFromChord),
               CurvePointNorm(Curve("cspline_ps"), hk_BlockStartFromChord),
               )

# PS

new_block_face(Point(p1.x, p1.y + yPerLowHGridBlockPitchStart * pitch - pitch, 0),
               Point(p2.x, p2.y + yPerLowHGridBlockPitchStart * pitch, 0),
               CurvePointNorm(Curve("cspline_ps"), hk_BlockStartFromChord),
               CurvePointNorm(Curve("cspline_ps"), vk_BlockStartFromChord)
               )

# SS

new_block_face(CurvePointNorm(Curve("cspline_ss"), hk_BlockStartFromChord),
               CurvePointNorm(Curve("cspline_ss"), vk_BlockStartFromChord),
               Point(p2.x, p2.y + yPerHighHGridBlockPitchStart * pitch, 0),
               Point(p1.x, p1.y + yPerHighHGridBlockPitchStart * pitch - pitch, 0)
               )

p3 = Point(CurvePointNorm(Curve("cspline_ss"), 0.5).x,
           CurvePointNorm(Curve("cspline_ss"), 0.5).y + yPerLowHGridBlockPitchStart * pitch, 0)

edge("Block_9", 1, 4).insert_vertex(0.5)
move_vertex(vertex("Block_9", 1, 4, 2), CurvePointNorm(Curve("cspline_ss"), 1))
edge("Block_8", 1, 1).insert_vertex(0.5)
move_vertex(vertex("Block_8", 1, 1, 2), p3)
edge("Block_12", 1, 2).insert_vertex(0.5)
move_vertex(vertex("Block_12", 1, 2, 2), p3)
edge("Block_10", 1, 2).insert_vertex(0.5)
move_vertex(vertex("Block_10", 1, 2, 2), CurvePointNorm(Curve("cspline_ss"), 0))
connect_whole_grid("ALL", 1E-06)

# =============================================================================
# Setze punktevertielungen
# =============================================================================


segment("Block_7", 1, 4, 1).set_number_of_points(
    int(factor * segment("Block_7", 1, 4, 1).get_discrete_length() / delta_i) + 1, 1)
segment("Block_5", 1, 4, 1).set_number_of_points(
    int(factor * segment("Block_5", 1, 4, 1).get_discrete_length() / delta_i) + 1, 1)
segment("Block_6", 1, 4, 1).set_number_of_points(
    int(factor * segment("Block_6", 1, 4, 1).get_discrete_length() / delta_i) + 1, 1)
segment("Block_2", 1, 3, 1).set_number_of_points(
    int(factor * segment("Block_2", 1, 3, 1).get_discrete_length() / delta_i) + 1, 1)
segment("Block_1", 1, 1, 1).set_number_of_points(
    int(factor * segment("Block_1", 1, 1, 1).get_discrete_length() / delta_i) + 1, 1)
segment("Block_4", 1, 1, 1).set_number_of_points(
    int(factor * segment("Block_4", 1, 1, 1).get_discrete_length() / delta_i) + 1, 1)
segment("Block_5", 1, 1, 1).set_number_of_points(
    int(factor * segment("Block_5", 1, 1, 1).get_discrete_length() / delta_i) + 1, 1)
segment("Block_8", 1, 2, 1).set_number_of_points(
    int(factor * segment("Block_8", 1, 2, 1).get_discrete_length() / delta_i) + 1, 1)
segment("Block_10", 1, 3, 1).set_number_of_points(
    int(factor * segment("Block_10", 1, 3, 1).get_discrete_length() / delta_i) + 1, 1)

segment("Block_1", 1, 1, 1).cluster_uniform()
segment("Block_4", 1, 1, 1).cluster_uniform()
segment("Block_5", 1, 1, 1).cluster_uniform()
segment("Block_5", 1, 4, 1).cluster_uniform()
segment("Block_6", 1, 4, 1).cluster_uniform()
segment("Block_7", 1, 4, 1).cluster_uniform()
segment("Block_7", 1, 2, 1).cluster_uniform()
segment("Block_8", 1, 2, 1).cluster_uniform()
segment("Block_3", 1, 2, 1).cluster_uniform()
segment("Block_3", 1, 3, 1).cluster_uniform()
segment("Block_2", 1, 3, 1).cluster_uniform()
segment("Block_1", 1, 3, 1).cluster_uniform()

segment("Block_10", 1, 2, 1).cluster_tanh(cellwidthcoeff * delta_i, cellwidthcoeff * delta_i)
segment("Block_9", 1, 4, 1).cluster_tanh(cellwidthcoeff * delta_i, cellwidthcoeff * delta_i)
segment("Block_12", 1, 1, 1).cluster_both_ends(cellwidthcoeff * delta_i)
segment("Block_11", 1, 2, 1).cluster_both_ends(cellwidthcoeff * delta_i)


# =============================================================================
# Glaettung
# =============================================================================

create_face_group("midspan")
face_group("midspan").add_face(face("Block_1", 1))
face_group("midspan").add_face(face("Block_2", 1))
face_group("midspan").add_face(face("Block_3", 1))
face_group("midspan").add_face(face("Block_4", 1))
face_group("midspan").add_face(face("Block_9", 1))
face_group("midspan").add_face(face("Block_11", 1))
face_group("midspan").add_face(face("Block_12", 1))
face_group("midspan").add_face(face("Block_8", 1))
face_group("midspan").add_face(face("Block_10", 1))
face_group("midspan").add_face(face("Block_5", 1))
face_group("midspan").add_face(face("Block_6", 1))
face_group("midspan").add_face(face("Block_7", 1))

create_segment_group("solid")
segment_group("solid").add_segment(segment("Block_11", 1, 2, 1))
segment_group("solid").add_segment(segment("Block_12", 1, 1, 1))
segment_group("solid").add_segment(segment("Block_9", 1, 4, 1))
segment_group("solid").add_segment(segment("Block_10", 1, 2, 1))
segment_group("solid").smoother_bc(0, first_cell_width, exp_ratio, layers)

create_segment_group("fixed")
segment_group("fixed").add_segment(segment("Block_1", 1, 1, 1))
segment_group("fixed").add_segment(segment("Block_4", 1, 1, 1))
segment_group("fixed").add_segment(segment("Block_5", 1, 1, 1))
segment_group("fixed").add_segment(segment("Block_5", 1, 4, 1))
segment_group("fixed").add_segment(segment("Block_6", 1, 4, 1))
segment_group("fixed").add_segment(segment("Block_7", 1, 4, 1))
segment_group("fixed").add_segment(segment("Block_7", 1, 2, 1))
segment_group("fixed").add_segment(segment("Block_8", 1, 2, 1))
segment_group("fixed").add_segment(segment("Block_3", 1, 2, 1))
segment_group("fixed").add_segment(segment("Block_3", 1, 3, 1))
segment_group("fixed").add_segment(segment("Block_2", 1, 3, 1))
segment_group("fixed").add_segment(segment("Block_1", 1, 3, 1))
segment_group("fixed").smoother_bc(3)


face_group("midspan").planar_smoothing(5000)

# =============================================================================
# Extrudieren
# =============================================================================


block_by_face_extrusion(face("Block_1", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_2", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_3", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_4", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_5", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_6", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_7", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_8", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_9", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_10", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_11", 1), Vector(0, 0, extrudeLength), 1, 1)
block_by_face_extrusion(face("Block_12", 1), Vector(0, 0, extrudeLength), 1, 1)

patch("Block_1", 3, 1).set_type("SOL")
patch("Block_1", 2, 1).set_type("SOL")
patch("Block_1", 5, 1).set_type("SOL")
patch("Block_2", 1, 1).set_type("SOL")
patch("Block_2", 2, 1).set_type("SOL")
patch("Block_2", 5, 1).set_type("SOL")
patch("Block_3", 1, 1).set_type("SOL")
patch("Block_3", 2, 1).set_type("SOL")
patch("Block_3", 4, 1).set_type("SOL")
patch("Block_3", 5, 1).set_type("SOL")
patch("Block_4", 1, 1).set_type("SOL")
patch("Block_4", 2, 1).set_type("SOL")
patch("Block_4", 3, 1).set_type("SOL")
patch("Block_5", 1, 1).set_type("SOL")
patch("Block_5", 2, 1).set_type("SOL")
patch("Block_5", 3, 1).set_type("SOL")
patch("Block_5", 6, 1).set_type("SOL")
patch("Block_6", 1, 1).set_type("SOL")
patch("Block_6", 2, 1).set_type("SOL")
patch("Block_6", 6, 1).set_type("SOL")
patch("Block_7", 1, 1).set_type("SOL")
patch("Block_7", 2, 1).set_type("SOL")
patch("Block_7", 4, 1).set_type("SOL")
patch("Block_7", 6, 1).set_type("SOL")
patch("Block_8", 1, 1).set_type("SOL")
patch("Block_8", 2, 1).set_type("SOL")
patch("Block_8", 4, 1).set_type("SOL")
patch("Block_9", 1, 1).set_type("SOL")
patch("Block_9", 2, 1).set_type("SOL")
patch("Block_9", 6, 1).set_type("SOL")
patch("Block_10", 1, 1).set_type("SOL")
patch("Block_10", 2, 1).set_type("SOL")
patch("Block_10", 4, 1).set_type("SOL")
patch("Block_11", 1, 1).set_type("SOL")
patch("Block_11", 2, 1).set_type("SOL")
patch("Block_11", 4, 1).set_type("SOL")
patch("Block_12", 1, 1).set_type("SOL")
patch("Block_12", 2, 1).set_type("SOL")
patch("Block_1", 1, 1).set_type("SOL")
patch("Block_12", 3, 1).set_type("SOL")

# search_connections(1E-007)


segment("Block_1", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_2", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_3", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_4", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_5", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_6", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_7", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_8", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_9", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_10", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_11", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)
segment("Block_12", 3, 3, 1).set_number_of_points(int(extrudeNodes * factor), 0)


search_connections(1E-007)

save_project(os.path.join(script_path, 'mesh.igg'))

export_FLUENT(os.path.join(script_path,"fluent.msh"))
