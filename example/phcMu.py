import sys
sys.path.insert(0, "../src")
import mx
import os
import math

baseName = "phcMu"
dim = 3
l = 1.0
o = 0.0
r = 16 #resolution in each direction

# set the grid
grid = mx.Grid(dim * [r], dim * [o], dim * [l])

# boundary conditions
bcs = mx.BoundaryConditions()
bcs.setLowerBCs(['periodic', 'periodic', 'periodic'])
bcs.setUpperBCs(['periodic', 'periodic', 'periodic'])
factor = 0.0e-8
bcs.setPhaseShifts([factor*math.pi/4., factor*math.pi/3., factor*math.pi])
#bcs.setPhaseShifts([0, 0, factor*math.pi])

sph = mx.Sphere(0.37, [l/2.,l/2.,l/2.])
mu = mx.Mu(sph, muDiag=[9.4, 9.4, 11.6], losstan=0.e-1)
#diel = mx.Dielectric(sph, epsDiag=[1.0, 1.0, 1.0])

# setup simulation
sim = mx.Simulation(baseName, dim)
sim.setGrid(grid)
sim.addMu(mu)
sim.setBCs(bcs)

# create new eigensolver object
eig = mx.Eigensolver()
eig.setParams({"nev": 10, "basis": 30})
#eig.setParams({"nev": 20, "basis": 30, "shift": -2.0})

lin = mx.LinearSolver()
lin.setParams({"sweeps": 2, "type": "gmres",
  "prec type": "amg",
  "smoother": "Chebyshev",
  "levels": 10, "basis": 20})

eig.setLinearSolver(lin)
sim.setSolver(eig)

sim.write()

os.system("../build/src/maxwell --infile=" + baseName + ".mx")
