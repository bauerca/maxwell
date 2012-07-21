import sys
sys.path.insert(0, "../src")
import mx
import os

baseName = "phc-sapph-r0.37"
dim = 3
l = 1.0
o = -0.5

# the object
sph = mx.Sphere(0.37, [0.0, 0.0, 0.0])
epsDiag = [10.225, 10.225, 9.95]
epsOffDiag = [0.67360967926537398, -0.67360967926537398, -0.825]
diel = mx.Dielectric(sph, epsDiag=epsDiag, epsOffDiag=epsOffDiag)

# the solver
solver = mx.Eigensolver()


#resolutions = range(8, 32)
resolutions = [36]
print resolutions

for r in resolutions:
  name = baseName + "-%.3d" % r

  grid = mx.Grid(dim * [r], dim * [o], dim * [l])
  sim = mx.Simulation(name, dim)

  sim.setGrid(grid)
  sim.addDielectric(diel)
  #sim.setLinearSolverParameters({"sweeps": 5, "type": "gmg-gmres", "levels": 3, "basis": 20})
  sim.setSolver(solver)
  sim.write()

  #os.system("mpirun -machinefile ./nodes -np 4 src/maxwell --infile=%s" % (name + ".mx"))
  #os.system("src/maxwell --infile=%s" % (name + ".mx"))
  #os.system("cp mxEigenfrequenciesReal.h5 eigfreqs-fit-phc-sapph/eigfreqs%.3d.h5" % r)

