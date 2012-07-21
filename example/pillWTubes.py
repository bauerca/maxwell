import sys
sys.path.insert(0, "../src")
import mx
import os
import math

baseName = "pillWTubes"
dim = 3

# units
mm = 1.e-3 # lengths in millimeters
ghz = 1.e9 # frequencies in GHz

# target frequency for TM_010 mode
freq010 = 12.0*ghz
omega010 = 2.*math.pi*freq010

# cavity iris dimensions
irisR = 3.15*mm
irisT = 1.67*mm

def getCavRadius(irisR, irisT, omega010, lzCav, syncPhAdv):
  """Use perturbation theory to predict the irises' effect on the TM_010 
  resonant frequency, and set the cavity radius to give the target
  TM_010 frequency.

  Parameters:
    irisR : iris radius
    irisT : iris thickness
    omega010 : target TM_010 angular frequency
    lzCav : length of pillbox part of cavity
    syncPhAdv : phase advance between cells for synchronous TM_010 mode at
                frequency omega010
  """
  
  # Pert theory says:
  #   a*omega_c^3 + omega_c - omega = 0
  #
  # where omega_c is the frequency of the pillbox cavity that will
  # give the target frequency, omega, once irises are added;
  #
  #   a = tau * \frac{irisR^3}{3 L c^2} * (1 - cos(phi)*e^(-alpha*irisT));
  #
  # where phi is the synchronous phase advance for the TM_010 mode, alpha is
  alpha = math.sqrt((2.405/irisR)**2 - (omega010/mx.lightspeed)**2)
  #
  #   tau = \frac{c^2 eps0 L E0^2}{omega_c^2 U}
  #       = \frac{c^2 E0^2}{omega_c^2 \int E^2 dx dy}
  #
  # where E0 is TM010 electric field on cavity axis, omega_c is closed
  # cavity frequency, and U is the closed cavity
  # energy.
  tau = 2. / (math.pi * (2.405 * 0.5191)**2)

  # solve the cubic equation:
  tmp1 = tau * (irisR**3 / (3.*lzCav*mx.lightspeed**2))
  tmp1 *= 1. - math.cos(syncPhAdv) * math.exp(-alpha * irisT)

  tmp2 = (math.sqrt(3.*27.*tmp1**4 * omega010**2 + 12.*tmp1**3) + 9. * tmp1**2 * omega010)**(1./3.)
  omegaPill = tmp2 / (2.**(1./3.) * 3.**(2./3.) * tmp1) - (2./3.)**(1./3.) / tmp2

  return mx.lightspeed * 2.405 / omegaPill

# synchronous phase advance (sets cavity length based on speed of light
# particle and TM_010 mode frequency)
syncPhAdv = 2.0*math.pi / 3.0
lz = mx.lightspeed * syncPhAdv / omega010
lzCav = lz - irisT

# phase advance (sets phase shift for periodic bcs in z-dir)
phAdv = 2.0 * math.pi / 3.0

# cavity radius
R = getCavRadius(irisR, irisT, omega010, lzCav, syncPhAdv)

# print cavity parameters
print "Periodic pillbox cavity parameters:"
print "  Iris radius = %.3f mm" % (irisR / mm)
print "  Iris thickness = %.3f mm" % (irisT / mm)
print "  Cavity radius = %.3f mm" % (R / mm)
print "  Cavity length = %.3f mm" % (lz / mm)
print "  Target TM_010 frequency = %.3f GHz" % (freq010 / ghz)
print "  Sync phase advance = %.3f Pi" % (syncPhAdv / math.pi)
print "  Simulated phase advance = %.3f Pi" % (phAdv / math.pi)

##########################
#  Geometry construction
##########################

# build the pillbox part of the cavity (no irises)
cyl = mx.Cylinder(R, mx.zhat, mx.zzz)
caps = mx.Slab(lzCav, mx.zhat, [0,0,0.5*lz])
pill = mx.ShapeIntersection([cyl, caps])

# build the irises (this part looks like a corrugated cylinder)
irisCyl = mx.Cylinder(irisR + 0.5*irisT, mx.zhat, mx.zzz)
torus1 = mx.Torus(irisR + 0.5*irisT, 0.5*irisT, mx.zhat, mx.zzz)
torus1.invert()
torus2 = mx.Torus(irisR + 0.5*irisT, 0.5*irisT, mx.zhat, [0,0,lz])
torus2.invert()
irisTube = mx.ShapeIntersection([irisCyl, torus1, torus2])

# combine
cav = mx.ShapeUnion([pill, irisTube])

##########################





##########################
# Grid setup
##########################

d = irisT / 4. # nominal cell size
lx = ly = 2.*R + 4.*d # 2-cell padding at transverse edge of cavity

nz = int(lz / d) + 1 
nx = ny = int(lx / d) + 1

origin = [-0.5*lx, -0.5*ly, 0.0]

# set the grid
grid = mx.Grid([nx, ny, nz], origin, [lx, ly, lz])

##########################




##########################
# Boundary conditions
##########################

bcs = mx.BoundaryConditions() # defaults to periodic in all directions
#bcs.setUpperBCs(["pec","pec","pec"])
#bcs.setLowerBCs(["pec","pec","pec"])
bcs.setPhaseShifts([0, 0, phAdv])
#bcs.setPhaseShifts([0, 0, 0])

##########################



##########################
# Simulation setup
##########################

sim = mx.Simulation(baseName, dim)
sim.setGrid(grid)
sim.addPEC(cav)
sim.setBCs(bcs)
sim.setEigensolverParameters({"nev": 10, "basis": 20})
sim.setLinearSolverParameters({"sweeps": 2, "type": "gmres", "prec type": "amg",
  "levels": 2, "basis": 20, "smoother": "Chebyshev"})
sim.write()

