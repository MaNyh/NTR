import sys
import os
import pickle

def read_pickle_args(path):
    print("reading yaml-dictionary ", path)
    filepath = os.path.join(path, "igg_args.pkl")
    with open(filepath, "rb") as Fobj:
        dict = pickle.load(Fobj)
    return dict


script_path = os.path.dirname(os.path.abspath(__file__))
package_path = os.path.join(script_path, "../../..", "..")
tmp_path = os.path.join(script_path)
args = read_pickle_args(tmp_path)

pointcloudfile = args["pointcloudfile"]
case_path = args["case_path"]
add_path = args["add_path"]

sys.path.append(add_path)
from NTR.utils.externals.tecplot.tecplot_functions import openTecplotFile



yPerLowHGridBlockPitchStart = args["yPerLowHGridBlockPitchStart"]
yPerHighHGridBlockPitchStart = args["yPerHighHGridBlockPitchStart"]
blockStartFromChord = args["blockStartFromChord"]
factor = args["factor"]
ogrid_factor = args["ogrid_factor"]
delta_i = args["delta_i"] / factor
cellwidthcoeff = args["cellwidthcoeff"]
streamline_nodedensity_factor = args["streamline_nodedensity_factor"]
first_cell_width = args["first_cell_width"] / factor

le_firstcellheight_coeff = args["le_firstcellheight_coeff"]
te_firstcellheight_coeff = args["te_firstcellheight_coeff"]

shift_vk_block_xaxiscoeff = args["shift_vk_block_xaxiscoeff"]
shift_hk_block_xaxiscoeff = args["shift_hk_block_xaxiscoeff"]

hk_ps_shift = args["hk_ps_shift"]
hk_ss_shift = args["hk_ss_shift"]
vk_ps_shift = args["vk_ps_shift"]
vk_ss_shift = args["vk_ss_shift"]

exp_ratio = args["exp_ratio"]

layers = args["layers"]
layers = int(layers * factor)

extrudeLength = args["extrudeLength"]
extrudeNodes = int(extrudeLength/delta_i*factor)#int(factor * args["extrudeNodes"])

smoothing_iterations = args["smoothing"]

save_project_path = args["save_project"]
save_fluent_path = args["save_fluent"]

print("factor: " + str(factor))
print("delta_i: " + str(delta_i))
print("first_cell_width: " + str(first_cell_width))
print("exp_ratio: " + str(exp_ratio))
print("layers: " + str(layers))
print("extrudeLength: " + str(extrudeLength))
print("extrudeNodes: " + str(extrudeNodes))

data = openTecplotFile(pointcloudfile)


def iggsplines_from_data():
    global pitch
    x_ss, y_ss = data[0][0], data[0][1]
    x_ps, y_ps = data[1][0], data[1][1]
    x_lower, y_lower = data[2][0], data[2][1]
    x_upper, y_upper = data[3][0], data[3][1]

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


blocks = ["Block_1", "Block_2", "Block_3", "Block_4",
          "Block_5", "Block_6", "Block_7", "Block_8",
          "Block_9", "Block_10", "Block_11", "Block_12"]


def extrude_to_3d():
    # =============================================================================
    # Extrudieren
    # =============================================================================
    for b in blocks:
        block_by_face_extrusion(face(b, 1), Vector(0, 0, extrudeLength), 1, 1)

    for b in blocks:
        segment(b, 3, 3, 1).set_number_of_points(extrudeNodes, 0)


def set_patches():
    patch("Block_1", 1, 1).set_type("SOL")
    patch("Block_1", 2, 1).set_type("SOL")
    patch("Block_1", 3, 1).set_type("SOL")
    patch("Block_1", 5, 1).set_type("SOL")

    patch("Block_2", 1, 1).set_type("SOL")
    patch("Block_2", 2, 1).set_type("SOL")
    patch("Block_2", 5, 1).set_type("SOL")

    patch("Block_3", 1, 1).set_type("SOL")
    patch("Block_3", 2, 1).set_type("SOL")
    patch("Block_3", 4, 1).set_type("SOL")
    patch("Block_3", 5, 1).set_type("SOL")

    patch("Block_4", 2, 1).set_type("SOL")
    patch("Block_4", 1, 1).set_type("SOL")
    patch("Block_4", 3, 1).set_type("SOL")

    patch("Block_5", 3, 1).set_type("SOL")
    patch("Block_5", 1, 1).set_type("SOL")
    patch("Block_5", 2, 1).set_type("SOL")
    patch("Block_5", 6, 1).set_type("SOL")

    patch("Block_6", 1, 1).set_type("SOL")
    patch("Block_6", 2, 1).set_type("SOL")
    patch("Block_6", 6, 1).set_type("SOL")

    patch("Block_7", 2, 1).set_type("SOL")
    patch("Block_7", 1, 1).set_type("SOL")
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
    patch("Block_12", 3, 1).set_type("SOL")
    connect_whole_grid("PATCHES", 1E-06)


def smooth_2d_mesh():
    # =============================================================================
    # Glaettung
    # =============================================================================
    create_face_group("midspan")
    for b in blocks:
        face_group("midspan").add_face(face(b, 1))

    create_segment_group("solid")
    segment_group("solid").add_segment(segment("Block_11", 1, 2, 1))
    segment_group("solid").add_segment(segment("Block_12", 1, 1, 1))
    segment_group("solid").smoother_bc(0, first_cell_width, exp_ratio, layers)


    create_segment_group("trailing")
    segment_group("trailing").add_segment(segment("Block_10", 1, 2, 1))
    segment_group("trailing").smoother_bc(0, first_cell_width*te_firstcellheight_coeff, exp_ratio, layers)

    create_segment_group("leading")
    segment_group("leading").add_segment(segment("Block_9", 1, 4, 1))

    segment_group("leading").smoother_bc(0, first_cell_width * le_firstcellheight_coeff, exp_ratio, layers)

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

    interval = 500
    smooth_runs = int(smoothing_iterations / interval)
    for i in range(smooth_runs):
        face_group("midspan").planar_smoothing(interval)


def set_nodedistribution():
    # =============================================================================
    # Setze punktevertielungen
    # =============================================================================

    segment("Block_7", 1, 4, 1).set_number_of_points(
        int(factor * segment("Block_7", 1, 4, 1).get_discrete_length() / delta_i * streamline_nodedensity_factor) + 1,
        1)
    segment("Block_5", 1, 4, 1).set_number_of_points(
        int(factor * segment("Block_5", 1, 4, 1).get_discrete_length() / delta_i * streamline_nodedensity_factor) + 1,
        1)
    segment("Block_6", 1, 4, 1).set_number_of_points(
        int(factor * segment("Block_6", 1, 4, 1).get_discrete_length() / delta_i * streamline_nodedensity_factor) + 1,
        1)
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
        int(ogrid_factor * factor * segment("Block_10", 1, 3, 1).get_discrete_length() / delta_i) + 1, 1)
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

    segment("Block_10", 1, 2, 1).cluster_tanh(cellwidthcoeff / factor**2, cellwidthcoeff / factor**2)
    segment("Block_9", 1, 4, 1).cluster_tanh(cellwidthcoeff / factor**2, cellwidthcoeff / factor**2)
    segment("Block_12", 1, 1,1).cluster_uniform()
    segment("Block_11", 1, 2,1).cluster_uniform()


def set_blocks():
    # =============================================================================
    # Blockerstellung
    # =============================================================================

    # VK-Blockgrenze
    p1 = CurvePointNorm(Curve("cspline_peri_upper"),
                        Curve("cspline_peri_upper").calc_normalize(
                            Curve("cspline_peri_upper").project_point(Point(CurvePointNorm(Curve("cspline_ps"), 0.0).x,
                                                                            CurvePointNorm(Curve("cspline_ps"),
                                                                                           0.0).y + shift_vk_block_xaxiscoeff * pitch,
                                                                            0))))

    # HK-Blockgrenze
    p2 = CurvePointNorm(Curve("cspline_peri_lower"),
                        Curve("cspline_peri_lower").calc_normalize(
                            Curve("cspline_peri_lower").project_point(Point(CurvePointNorm(Curve("cspline_ps"), 1.0).x,
                                                                            CurvePointNorm(Curve("cspline_ps"),
                                                                                           0.0).y + shift_hk_block_xaxiscoeff * pitch,
                                                                            1.0))) + 0.05)
    p3 = Point(CurvePointNorm(Curve("cspline_ss"), 0.5).x,
               CurvePointNorm(Curve("cspline_ss"), 0.5).y + yPerLowHGridBlockPitchStart * pitch, 0)

    #VK
    pt1 = Point(p1.x, p1.y + yPerHighHGridBlockPitchStart * pitch - pitch, 0)  # vk_ss
    pt2 = CurvePointNorm(Curve("cspline_ss"), (1-blockStartFromChord)*vk_ss_shift)  # vk_ssBOUNDARY
    pt3 = CurvePointNorm(Curve("cspline_ps"), blockStartFromChord*vk_ps_shift)  # vk_psBOUNDARY
    pt4 = Point(p1.x, p1.y + yPerLowHGridBlockPitchStart * pitch - pitch, 0) # vk_ss
    #HK
    pt5 = Point(p2.x, p2.y + yPerLowHGridBlockPitchStart * pitch, 0)
    pt6 = Point(p2.x, p2.y + yPerHighHGridBlockPitchStart * pitch, 0)
    pt7 = CurvePointNorm(Curve("cspline_ss"), blockStartFromChord*hk_ss_shift )
    pt8 = CurvePointNorm(Curve("cspline_ps"), (1-blockStartFromChord)*hk_ps_shift)

    pt9 = CurvePointNorm(Curve("cspline_peri_lower"), 0.0)
    pt10 = CurvePointNorm(Curve("cspline_peri_lower"),Curve("cspline_peri_lower").calc_normalize(Curve("cspline_peri_lower").project_point(Point(p1.x, p1.y - pitch, 0))))
    pt11 = CurvePointNorm(Curve("cspline_inlet"), yPerLowHGridBlockPitchStart)

    pt12 = CurvePointNorm(Curve("cspline_inlet"), yPerHighHGridBlockPitchStart)

    pt13 = CurvePointNorm(Curve("cspline_peri_upper"), 0.0)

    pt14 = CurvePointNorm(Curve("cspline_peri_lower"), Curve("cspline_peri_lower").calc_normalize(Curve("cspline_peri_lower").project_point(Point(p1.x, p1.y - pitch, 0))))
    pt15 = CurvePointNorm(Curve("cspline_peri_lower"), Curve("cspline_peri_lower").calc_normalize(Curve("cspline_peri_lower").project_point(Point(p2.x, p2.y, 0))))

    pt16 = CurvePointNorm(Curve("cspline_peri_lower"), Curve("cspline_peri_lower").calc_normalize(Curve("cspline_peri_lower").project_point(Point(p2.x, p2.y, 0))))
    pt17 = CurvePointNorm(Curve("cspline_peri_lower"), 1.0)
    pt18 = CurvePointNorm(Curve("cspline_outlet"), yPerLowHGridBlockPitchStart)
    pt19 = CurvePointNorm(Curve("cspline_outlet"), yPerHighHGridBlockPitchStart)
    pt20 = CurvePointNorm(Curve("cspline_peri_upper"), 1.0)
    pt21 = CurvePointNorm(Curve("cspline_peri_upper"), Curve("cspline_peri_upper").calc_normalize(Curve("cspline_peri_upper").project_point(Point(p2.x, p2.y + pitch, 0))))

    """
    print(pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11,
          pt12, pt13, pt14, pt15, pt16, pt17, pt18, pt19, pt20, pt21, p1, p2, p3 )
    """

    new_block_face(pt9, pt10, pt4, pt11)
    new_block_face(pt11, pt4, pt1, pt12)
    new_block_face(pt12, pt1, p1, pt13)
    new_block_face(pt14, pt15, pt5, pt4)
    new_block_face(pt16, pt17, pt18, pt5)
    new_block_face(pt5,pt18, pt19, pt6)
    new_block_face(pt6, pt19, pt20, pt21)
    new_block_face(pt1, pt6, pt21, p1)

    # O-Netz

    # Vk
    new_block_face(pt4, pt3, pt2, pt1)

    # Hk
    new_block_face(pt5, pt6, pt7, pt8)

    # PS
    new_block_face(pt4, pt5, pt8, pt3)

    # SS
    new_block_face(pt2, pt7, pt6, pt1)

    edge("Block_9", 1, 4).insert_vertex(0.5)
    move_vertex(vertex("Block_9", 1, 4, 2), CurvePointNorm(Curve("cspline_ss"), 1))
    edge("Block_8", 1, 1).insert_vertex(0.5)
    move_vertex(vertex("Block_8", 1, 1, 2), p3)
    edge("Block_12", 1, 2).insert_vertex(0.5)
    move_vertex(vertex("Block_12", 1, 2, 2), p3)
    edge("Block_10", 1, 2).insert_vertex(0.5)
    move_vertex(vertex("Block_10", 1, 2, 2), CurvePointNorm(Curve("cspline_ss"), 0))
    connect_whole_grid("ALL", 1E-06)


iggsplines_from_data()
set_blocks()
set_nodedistribution()
smooth_2d_mesh()
extrude_to_3d()
set_patches()

os.chdir(case_path)
save_project("mesh.igg")
export_FLUENT("fluent.msh")
