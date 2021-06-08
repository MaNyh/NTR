import mFunc

import pyvista as pv

import numpy as np
from scipy.spatial.distance import pdist, squareform



# files = ["2D_Profile_c.dat","2D_Profile_m.dat","Q3D_Profile_c.dat","Q3D_Profile_m.dat"]
Blade = r'..\examples\cascade_postprocessing\gwk_verdichtercascade_les\VTK\BLADE\BLADE_81000.vtk'
Blade_Surf = pv.PolyData(Blade)
Blade_Surf.clear_arrays()
Blade_PTS = Blade_Surf.slice(normal="z", origin=(0, 0, 0.01))

camberLineResolution = 300
resolutionCutLineMiddle = 300
scanRangeFromCamberLength = 0.8
scaleLineSec = camberLineResolution / 10


def linesFromPoints(points):
    poly = pv.PolyData()
    poly.points = points
    cells = np.full((len(points) - 1, 3), 2, dtype=np.int_)
    cells[:, 1] = np.arange(0, len(points) - 1, dtype=np.int_)
    cells[:, 2] = np.arange(1, len(points), dtype=np.int_)
    poly.lines = cells
    return poly


def findMiddle(input_list):
    middle = float(len(input_list)) / 2
    if middle % 2 != 0:
        return input_list[int(middle - .5)]
    else:
        p0 = input_list[int(middle) - 1]
        p1 = input_list[int(middle)]
        return p0 + 0.5 * (p1 - p0)


def distance(pt1, pt2):
    diff = pt1 - pt2
    return diff[0] ** 2 + diff[1] ** 2 + diff[2] ** 2


def lineToVector(pvLine):
    start = pvLine.points[0]
    end = pvLine.points[-1]

    return end - start


# sehnenlinie bestimmen
# res = max(combinations(fPts, 2), key = lambda sub: distance(sub[0],sub[1]))

A = Blade_PTS.points
D = squareform(pdist(A))
N = np.max(D)
I = np.argmax(D)
I_row, I_col = np.unravel_index(I, D.shape)

pt1 = A[I_row]
pt2 = A[I_col]

allpts = list(Blade_PTS.points)

c = pv.Line(pt1, pt2, camberLineResolution)
c.points[::, 2] = 0

plotter = pv.Plotter()

# plotter = pvq.BackgroundPlotter()
poly = pv.PolyData(allpts)
# oly.plot()
face = Blade_Surf  # poly.delaunay_2d()
surface = face.extract_feature_edges()
surface.points[::, 2] = 0
plotter.add_mesh(face, point_size=30, color="b")
#####################################plotter.add_mesh(c)
# plotter.add_mesh(Blade_PTS)


centerLine = []
for idx, pt in enumerate(c.points):
    # erzeuge "cutline", diese soll sich mit bladeSS und bladePS schneiden! (surface)
    cutLine = pv.Line(scanRangeFromCamberLength * (pt1 - pt2), (0, 0, 0), resolutionCutLineMiddle)
    cutLine.translate(-cutLine.points[0] / 2)
    cutLine.rotate_z(90)
    cutLine.translate(pt)

    intersect = surface.slice_along_line(cutLine)
    if intersect.number_of_points > 1:
        centerLine.append(intersect.center)

beginDir = np.array(centerLine[1]) - np.array(centerLine[0])
endDir = np.array(centerLine[-2]) - np.array(centerLine[-1])

beginExt = -beginDir * 2 + centerLine[0]
endExt = -endDir * 2 + centerLine[-1]

centerLine = [beginExt] + centerLine + [endExt]
skeletalLine = linesFromPoints(np.array(centerLine))

######################
correctedChordStart = pv.Line(beginExt, centerLine[1], resolutionCutLineMiddle)
correctedChordEnd = pv.Line(centerLine[-2], endExt, resolutionCutLineMiddle)
interSectStart = surface.slice_along_line(correctedChordStart)
interSectEnd = surface.slice_along_line(correctedChordEnd)

correctedChordLine = pv.Line(interSectStart.points[-1], interSectEnd.points[0])

######################
# ermittelte skelettlinie
plotter.add_mesh(skeletalLine, color='w')
# neue chord
# plotter.add_mesh(correctedChordLine)

###GRÜN = TRAILINGEDGE
plotter.add_mesh(interSectStart, color='y', point_size=5)
###ROT = LEADINGEDGE
plotter.add_mesh(interSectEnd, color='r', point_size=5)

ssPts = []
psPts = []

for idx in range(skeletalLine.number_of_cells):
    workLine = skeletalLine.extract_cells(idx)
    p0 = workLine.points[0].copy()

    workLine.points -= p0
    workLine.points *= scaleLineSec
    workLine.rotate_z(90)
    workLine.points -= workLine.points[-1] * .5
    # workLine.points *= 10
    workLine.points += p0
    pts = workLine.points
    intersect = surface.slice_along_line(pv.Line(pts[1], pts[0], 2))
    if intersect.number_of_points > 0:
        # plotter.add_mesh(intersect)
        # plotter.add_mesh(workLine)
        ssPts.append(intersect.points[0])
        psPts.append(intersect.points[-1])

# pkt ps and ss blade
ssLine = pv.PolyData(np.array(ssPts))
psLine = pv.PolyData(np.array(psPts))

# geschwindigkeit absolut
# Ma = 0.25
# rho =
v = 68
u = np.array([v, 0, 0])
# anstellwinkel
alpha = 49.02

print("Winkel zwischen x-Achse und sehne")
chordXBeta = mFunc.angle_between(np.array([1, 0, 0]), lineToVector(c)) / (2 * np.pi) * 360
# LE:
scelXBeta2 = mFunc.angle_between(lineToVector(interSectStart), np.array([1, 0, 0])) / (2 * np.pi) * 360
# TE:
scelXBeta1 = mFunc.angle_between(lineToVector(interSectEnd), np.array([1, 0, 0])) / (2 * np.pi) * 360  # TE

print(chordXBeta)

xLine = pv.Line((0, 0, 0), (1, 0, 0))
xLine.points *= c.length / xLine.length
metal = pv.Line((0, 0, 0), interSectStart.points[1] - interSectStart.points[0])
metal.points *= c.length / metal.length
plotter.add_mesh(metal, color="r")
plotter.add_mesh(xLine, color="w")
uLine = metal.copy()
uLine.rotate_z(-chordXBeta + alpha)
u = uLine.points[1] - uLine.points[0]
uv = mFunc.vecDir(u) * v
print(uv)
plotter.add_mesh(uLine, color="b")
plotter.show()
