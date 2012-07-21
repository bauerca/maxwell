import os
import math
import sys

sys.path.insert(0, "../src")
import mx

baseName = "crabcav"
dim = 3

# crab cav parameters
#   user set:
numCells = 4
cellLen = 2.0 * 0.0192
cavRad = 0.04719
irisRad = 0.015
cavRho = 0.0136
irisRho = 0.00331

#   calculated:
irisTorMajRad = irisRad + irisRho
irisTorMinRad = irisRho
cavTorMajRad = cavRad - cavRho
cavTorMinRad = cavRho
rhoSum = cavRho + irisRho
rhoSum2 = rhoSum * rhoSum
radDiff = cavRad - irisRad
halfCellLen2 = 0.25 * cellLen * cellLen
diff2 = (radDiff - rhoSum)**2

cosTheta = (rhoSum - radDiff) * rhoSum
cosTheta += math.sqrt(halfCellLen2 * (diff2 - rhoSum2 + halfCellLen2))
cosTheta /= halfCellLen2 + diff2

theta = math.acos(cosTheta)
sinTheta = math.sqrt(1 - cosTheta * cosTheta)
tanTheta = sinTheta / cosTheta
cotTheta = 1.0 / tanTheta

coneVertOffset = 0.5 * cellLen - irisRho * sinTheta + (irisRad + irisRho * (1.0 - cosTheta)) * cotTheta



# build cavity using fundamental shapes
zhat = [0, 0, 1]
origin = [0, 0, 0]

irisTube = mx.Cylinder(irisRad + irisRho * (1.0 - cosTheta), zhat, origin)
irisTorus = mx.Torus(irisTorMajRad, irisTorMinRad, zhat, [0, 0, 0.5 * cellLen])

corrIrisTube = mx.ShapeSubtract(base = irisTube, subtract = irisTorus)

cavTube = mx.Cylinder(cavRad - cavRho * (1.0 - cosTheta), zhat, origin)
cavCone = mx.Cone(theta, zhat, [0, 0, coneVertOffset])
cavTorus = mx.Torus(cavTorMajRad, cavTorMinRad, zhat, origin)

preCav = mx.ShapeIntersection([cavCone, cavTube])
halfCell = mx.ShapeUnion([preCav, cavTorus, corrIrisTube])

fullCell = mx.ShapeMirror(halfCell, zhat, origin)

infCells = mx.ShapeRepeat(fullCell, origin, zhat, cellLen, numCells / 2, numCells / 2)

caps = mx.Slab(float(numCells) * cellLen, zhat, origin)

cav = mx.ShapeIntersection([caps, infCells])


# grid setup
cellRes = 10
pad = 2

delta = cellLen / float(cellRes)
nz = numCells * cellRes + 2 * pad
lz = float(nz) * delta
oz = -0.5 * lz

ny = nx = 2 * (int(math.ceil(cavRad / delta)) + pad)
ly = lx = float(nx) * delta
oy = ox = -0.5 * lx

grid = mx.Grid([nx, ny, nz], [ox, oy, oz], [lx, ly, lz])


# simulation setup

sim = mx.Simulation("crabcav.mx", 3)
sim.setGrid(grid)
sim.addPEC(cav)
sim.setLinearSolverParameters({"sweeps": 3})
sim.setEigensolverParameters({"nev": 15})
sim.write()

sys.exit()


resolutions = range(36, 40)
print resolutions

for r in resolutions:
  name = baseName + "-%.3d" % r

  grid = mx.Grid(dim * [r], dim * [o], dim * [l])
  sim = mx.Simulation(name + ".mx", dim)

  sim.setGrid(grid)
  sim.addDielectric(diel)
  sim.setLinearSolverParameters({"sweeps": 3})
  sim.write()

  os.system("mpirun -machinefile ./nodes -np 4 src/maxwell --infile=%s" % (name + ".mx"))
  os.system("cp mxEigenfrequenciesReal.h5 eigfreqs-fit-phc-sapph/eigfreqs%.3d.h5" % r)

