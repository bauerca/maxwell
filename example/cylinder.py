#!/usr/bin/env python
"""
@file    cylinder.py

@brief   To find the propagation modes for a cylindrical fiber

@version $Id$

Copyright &copy; 2017-2017, Tech-X Corporation, Boulder, CO.
"""

# Most dependent to least
import scipy
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

# Should be able to get scipy.constants, but since not
speed_of_light = 299792458.0

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

# Basic parameters from Zhu, Brown,
# "Full-vectorial finite-difference analysis of microstructured optical fibers"
# The calculation window is chosen to be the first quadrant of the fiber
# cross section with a computation window size of 6um by 6 um.
  fiberRadius = 3.e-6
  refractionIndex = 1.45
  n_eff = 1.438604 # Computed value
  wavelen = 1.5e-6
  relPermittivity = refractionIndex**2
  vacWavelen = wavelen*refractionIndex
  freq = speed_of_light/(wavelen/n_eff)
  print("Vacuum wavelength = %g, frequency = %g"%(vacWavelen, freq))

# Grid parameters
  resolution = 0.05 # Number of cells per wavelength
  endyz = 6.e-6
  bgnyz = -endyz
  lenyz = endyz - bgnyz
  dl = resolution*fiberRadius
  yzcells_half = int(endyz/dl)
  yzcells = 2*yzcells_half
  xcells = 2
  lenx = xcells*dl
  print("%d cells in the x direction, %d cells in the y and z directions."%
      (xcells, yzcells))
# set the grid
  dim = 3
  grid = mx.Grid([xcells, yzcells, yzcells], [0., bgnyz, bgnyz],
      [lenx, lenyz, lenyz])

# boundary conditions
  bcs = mx.BoundaryConditions()
  bcs.setLowerBCs(['periodic', 'pec', 'pec'])
  bcs.setUpperBCs(['periodic', 'pec', 'pec'])
  kay = n_eff*2.*math.pi/vacWavelen
  phaseShift = kay*lenx
  bcs.setPhaseShifts([phaseShift, 0., 0.])
  print("phase shift = %g."%(phaseShift))

  cylLen = lenx + 4.*dl
  cylStart = - 2.*dl
  cyl = mx.Cylinder(fiberRadius, [1., 0., 0.], [cylStart, 0., 0.])
# mu = mx.Mu(sph, muDiag=[9.4, 9.4, 11.6], losstan=0.e-1)
  diel = mx.Dielectric(cyl, epsDiag=[relPermittivity, relPermittivity,
      relPermittivity])

# setup simulation
  baseName = "cylinder" # Should come from args
  sim = mx.Simulation(baseName, dim)
  sim.setGrid(grid)
  sim.addDielectric(diel)
  sim.setBCs(bcs)

# create new eigensolver object
  eig = mx.Eigensolver()
  eig.setParams({"nev": 10, "basis": 30})
  prec = mx.AMG()
  lin = mx.LinearSolver(prec = prec)
# lin.setParams({"sweeps": 2, "type": "gmres",
  # "prec type": "amg",
  # "smoother": "Chebyshev",
  # "levels": 10, "basis": 20})
  eig.setLinearSolver(lin)
  sim.setSolver(eig)

# Write input file and run simulation
  sim.write()
  mxwl = "../../builds/maxwell/ser/src/maxwell"
  exline = mxwl + " --infile=" + baseName + ".mx"
  print(exline)
  if os.path.isfile(mxwl):
    os.system(exline)

if __name__ == '__main__':
  print("cylinder.py started.")
  main()

