#!/usr/bin/env python
"""
@file    cylinder.py

@brief   To find the propagation modes for a cylindrical fiber

@version $Id$

Copyright &copy; 2017-2017, Tech-X Corporation, Boulder, CO.
"""

# Most dependent to least
import math
import argparse
import sys
import os

# Add path to mx, etc
myexec = sys.argv[0]
mydir = os.path.dirname(os.path.realpath(myexec))
mytopdir = os.path.dirname(mydir)
# print("mytopdir = %s."%mytopdir)
for dr in [os.path.join(mytopdir, "share"), os.path.join(mytopdir, "bin"),
    os.path.join(mytopdir, "src")]:
  if os.path.exists(dr):
    sys.path.insert(0, dr)
# print("sys.path:")
# print(sys.path)
import mx

"""
baseName = "cylinder"
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
# mu = mx.Mu(sph, muDiag=[9.4, 9.4, 11.6], losstan=0.e-1)
diel = mx.Dielectric(sph, epsDiag=[2.25, 2.25, 2.25])

# setup simulation
sim = mx.Simulation(baseName, dim)
sim.setGrid(grid)
sim.addDielectric(diel)
sim.setBCs(bcs)

# create new eigensolver object
eig = mx.Eigensolver()
eig.setParams({"nev": 10, "basis": 30})
#eig.setParams({"nev": 20, "basis": 30, "shift": -2.0})

lin = mx.LinearSolver()

# lin.setParams({"sweeps": 2, "type": "gmres",
  # "prec type": "amg",
  # "smoother": "Chebyshev",
  # "levels": 10, "basis": 20})


eig.setLinearSolver(lin)
sim.setSolver(eig)

sim.write()
mxwl = "../../builds/maxwell/ser/src/maxwell"
exline = mxwl + " --infile=" + baseName + ".mx"
print(exline)
if os.path.isfile(mxwl):
  os.system(exline)
"""

#
# Main
#
def main():

# create the cmd line parser
  parser = argparse.ArgumentParser()

# Just examples for now
  parser.add_argument("-m", "--maximum-time",
      help="Maximum time (seconds) the job will run",
      action="store", dest="maxtime", type=int)
  parser.add_argument("-t", "--testing", help="test signal catching",
      action="store_true")
  (args, unkargs) = parser.parse_known_args()

if __name__ == '__main__':
  print("cylinder.py started.")
  main()
