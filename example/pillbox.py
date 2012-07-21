import sys
sys.path.insert(0, "../src")
import mx
import os

baseName = "pillbox"
pol = "TE"
dim = 3
l = 1.0
o = -0.5
r = 20 #resolution in each direction

# build the cavity
cyl = mx.Cylinder(0.4, [0,0,1], [0,0,0])
caps = mx.Slab(0.8, [0,0,1], [0,0,0])
cav = mx.ShapeIntersection([cyl, caps])
pec = mx.PEC(cav)

# set the grid
grid = mx.Grid(dim * [r], dim * [o], dim * [l])

##########################
# create new eigensolver object
##########################
eig = mx.Eigensolver()
eig.setParams({"nev": 2, "basis": 100, "spectrum": "LM"})
#eig.setParams({"nev": 100, "basis": 150, "shift": -1000})

lin = mx.LinearSolver(solveType="gmres", basis=40)

#prec = mx.AMG(sweeps=2)
prec = mx.ILU(dropTol=1.e-3, fill=10.)
lin.setPreconditioner(prec)

eig.setLinearSolver(lin)
##########################



##########################
# Simulation setup
##########################

sim = mx.Simulation(baseName, dim, pol)
sim.setGrid(grid)
sim.setSolver(eig)
sim.addPEC(pec)
sim.write()

os.system("../cmakebuild/src/maxwell --infile=" + baseName + ".mx")
#os.system("../build/src/maxwell --infile=" + baseName + ".mx")

