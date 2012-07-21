import sys
sys.path.insert(0, "../src")
import mx
import os
import math

baseName = "vacuum"
dim = 3
pol = "TM"
l = 1.0
o = 0.0
r = 16 #resolution in each direction

# set the grid
grid = mx.Grid(dim * [r], dim * [o], dim * [l])

# boundary conditions
bcs = mx.BoundaryConditions()
bcs.setLowerBCs(['periodic', 'periodic', 'periodic'])
bcs.setUpperBCs(['periodic', 'periodic', 'periodic'])
factor = 0.0e-2
bcs.setPhaseShifts([factor*math.pi/4., factor*math.pi/3., factor*math.pi])
#bcs.setPhaseShifts([0, 0, factor*math.pi])


##########################
# create new eigensolver object
##########################
eig = mx.Eigensolver()
#eig.setParams({"nev": 10, "basis": 100, "spectrum": "LM", "shift": 500.})
#eig.setParams({"nev": 100, "basis": 150, "shift": -1000})

#lin = mx.LinearSolver(solveType="gmres", basis=40)

#prec = mx.AMG(sweeps=2)
#prec = mx.ILU(dropTol=1.e-3, fill=10.)
#lin.setPreconditioner(prec)

#eig.setLinearSolver(lin)
##########################

##########################
# Simulation setup
##########################

sim = mx.Simulation(baseName, dim, pol)
sim.setGrid(grid)
sim.setBCs(bcs)
sim.setSolver(eig)
sim.write()

#os.system("../cmakebuild/src/maxwell --infile=" + baseName + ".mx")
#os.system("../build/src/maxwell --infile=" + baseName + ".mx")
