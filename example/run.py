import mx
import os
from sys import exit
import numpy
import math
#import dsph_in_msph as dm
import dsphmsph as dm
import tables

simName = "dsphmsph"
dim = 3
simPol = "TE"

# grid info
lx = ly = lz = 0.5
origin = 3 * [0.0]

# setup dielectric sphere
epsRad = 0.37
eps = 10.0
epsDiag = 3 * [eps]
epsOffDiag = 3 * [0.0]
epsSph = mx.Sphere(epsRad, origin)
diel = mx.Dielectric(epsDiag, epsOffDiag, epsSph)

# setup pec sphere
pecRad = 0.49
pecSph = mx.Sphere(pecRad, origin)

# eigensolver params
nev = 3
basis = 20
#eigparams = {"nev" : nev, "basis" : basis, "shift" : 0.5}
globEigParams = {"basis" : basis, "shift" : 0.5}

# global variables for field analysis
shellPad = 3 # distance from dielectric boundary to evaluate field errors (units of dx)
#shellPads = range(-6, -2) + range(3, 5)
#shellPads = range(3, 7)
shellPads = [-3, 3]

subdir = "energy_norm_min_res_32"



def enumBCs(options, nslots):
  nopts = len(options)
  N = int(math.pow(nopts, nslots))
  res = [['' for j in range(nslots)] for i in range(N)]
  print res
  for i in range(N): # combination index
    for j in range(nslots): # slot index for current combination
      res[i][j] = options[(i / int(math.pow(nopts, j))) % nopts]
  return res

#print enumBCs(["pmc", "pec"], 3)



#simModes(["pmc", "pmc", "pec"])
def getDegen(modeNames, modeFreqs):
  degens = []
  degen = []
  for i in range(1, len(modeFreqs)):
    match = (modeFreqs[i] == modeFreqs[i - 1])
    if match:
      if degen == []:
        degen.append(i - 1)
      degen.append(i)
      if i == len(modeFreqs) - 1:
        degens.append(degen)
    elif degen:
      degens.append(degen)
      degen = []
  return degens

#print getDegen([1, 2, 3, 4, 5], [1, 1, 3, 3, 3])
#exit()
  

def getShellPoints(n, minRes, thetaRes, pad):
  """
  Get a set of points on a shell just outside of inner dielectric sphere
    n:           simulation resolution (n = nx = ny = nz assumed)
    thetaRes:    resolution in the theta direction from 0 to pi/2 (dPhi = dTheta)
    numCellPad:  distance from dielectric interface in units of the grid cell x-length
  """
  n = float(n)
  dx = lx / n
  dTheta = 0.5 * math.pi / float(thetaRes)
  shellRad = epsRad + float(pad) * dx
  maxDx = (1. + 1.e-8) * lx / float(minRes) 
  thetaMin = math.asin(maxDx / shellRad / math.sqrt(2.0))
  thetaMax = math.acos(0.5 * maxDx / shellRad)
  points = []
  for theta in numpy.arange(thetaMin, thetaMax, dTheta):
    sinTh = math.sin(theta)
    phiMin = math.asin(0.5 * maxDx / (shellRad * sinTh))
    phiMax = math.acos(0.5 * maxDx / (shellRad * sinTh))
    for phi in numpy.arange(phiMin, phiMax, dTheta / sinTh):
      points.append([shellRad * math.sin(theta) * math.cos(phi),
                     shellRad * math.sin(theta) * math.sin(phi),
                     shellRad * math.cos(theta)])
  return points


def getVolumePoints(minRes, rRes, region):
  """
  Get a set of points for calculating volume norms. Points are restricted to the following
  sub-volumes: first, points must lie a distance 0.5dx away from all simulation domain
  boundaries, so that an interpolation doesn't try to use values outside the domain (Yee
  setup making life difficult again); second, points must lie far enough away from material
  boundaries, so that complicated interpolations are not required.
  """

  # when every resolution has the same bndry buffer
  maxDx = (1. + 1.e-8) * lx / float(minRes)
  dr = pecRad / float(rRes)

  # shell distances inside dielectric
  rmin = 0.5 * math.sqrt(3.0) * maxDx
  rmax = epsRad - 3.0 * maxDx
  rIn = numpy.arange(rmin, rmax, dr)

  # shell distances outside dielectric
  rmin = epsRad + 3.0 * maxDx
  rmax = pecRad - 3.0 * maxDx
  rOut = numpy.arange(rmin, rmax, dr)

  if region == "in":
    rs = rIn
  elif region == "out":
    rs = rOut
  else:
    rs = numpy.concatenate([rIn, rOut])

  points = []
  for r in rs:
    dTheta = math.acos(1.0 - 0.5 * (dr / r)**2)
    thetaMin = math.asin(maxDx / r / math.sqrt(2.0))
    thetaMax = math.acos(0.5 * maxDx / r)
    for theta in numpy.arange(thetaMin, thetaMax, dTheta):
      sinTh = math.sin(theta)
      dPhi = dTheta / sinTh
      phiMin = math.asin(0.5 * maxDx / (r * sinTh))
      phiMax = math.acos(0.5 * maxDx / (r * sinTh))
      for phi in numpy.arange(phiMin, phiMax, dPhi):
        points.append([r * math.sin(theta) * math.cos(phi),
                       r * math.sin(theta) * math.sin(phi),
                       r * math.cos(theta)])
  return points
        
  

def savePoints(points, dir):
  h5points = tables.openFile(dir + "shell_points.h5", 'w')
  h5points.createArray("/", "data", points)
  h5points.close()

def writeSim(n, lbcs, pointSets):
  sim = mx.Simulation(simName, dim)

  grid = mx.Grid(3 * [n], origin, [lx, ly, lz])

  output = mx.Output()
  for name, points in pointSets.iteritems():
    output.fieldValues(name, points)

  bcs = mx.BoundaryConditions()
  bcs.setUpperBCs(["pec", "pec", "pec"])
  bcs.setLowerBCs(lbcs)

  eigParams = globEigParams
  eigParams["nev"] = len(dm.octModes(lbcs, nev)[0])

  sim.setGrid(grid)
  sim.addDielectric(diel)
  sim.addPEC(pecSph)
  sim.setBCs(bcs)
  sim.setOutput(output)
  sim.setEigensolverParameters(eigParams)
  sim.write()


def simDir(resolution, lbcs):
  resDir = "./results_%.4d/" % resolution
  dir = "./results_%.4d/%s_%s_%s/" % (resolution, lbcs[0], lbcs[1], lbcs[2])
  if not os.path.exists(resDir):
    os.mkdir(resDir)
  if not os.path.exists(dir):
    os.mkdir(dir)
  return dir

def filesContaining(dir, substrs):
  files = os.listdir(dir)
  for s in substrs:
    files = [file for file in files if s in file]
  return files

def saveResults(n, lbcs):
  """
  Moves all mx output files into the new directory: './results_%.4d/' % n. Also, eigenmode h5 files 
  (and h5 files representing eigenmode values at selected points) are renamed using spherical multipole
  notation. For example, if 'vec00' refers to the TM010re mode (based on eigenvalue ordering and
  boundary conditions), the file,
    mxYeeElecField_eigenvectors_vec00_comp0.h5
  will be mv'ed to,
    ./results_####/mxYeeElecField_eigenvectors_TM010re_comp0.h5
  (Warning: this may get the labeling of degenerate modes incorrect. This will be taken into account
  during the analysis of the fields.)
    
    n = resolution
    lbcs = -x,-y,-z boundary conditions
  """

  dir = simDir(n, lbcs)
  #if dim == 3:
  #  mx.h5tovtk("./", "~/cerberus-mainline/h5utils/bin/")

  # get mode names and analytic frequencies that should be found for these boundary conditions
  # (sorted by increasing mode frequency)
  modeNames, modeFreqs = dm.octModes(lbcs, nev)

  # get field h5 output from mx
  fieldNames = filesContaining("./", ["Field", ".h5"])
  fieldNames.sort()

  # replace mx naming scheme with actual mode names (e.g. substring 'vec00' --> 'TM010re')
  for name in fieldNames:
    vecStr = name.split("_")[-1].rstrip(".h5")
    vecNo = int(vecStr.lstrip("vec"))
    newName = name.replace(vecStr, modeNames[vecNo])
    os.system("mv " + name + " " + dir + newName)

  # add mx frequencies to a list of all frequencies found at this resolution
  h5freq = tables.openFile("mxEigenfrequenciesReal.h5", 'r')
  freqs = h5freq.root.data.read()
  h5freq.close()
  freqDiff = numpy.abs(freqs - modeFreqs) / modeFreqs
  numpy.savetxt(dir + "relFreqDiffs.txt", freqDiff)
  os.system("mv mxEigenfrequencies* " + dir)

  # save true frequencies
  f = open(dir + "trueFreqs.txt", 'w')
  for modeName, modeFreq in zip(modeNames, modeFreqs):
    f.write(modeName + ", %.8e\n" % modeFreq)
  f.close()

def mxFieldValuesAtPoints(resolution, lbcs, type, modename, pointSet, pointSetName):
  dir = simDir(resolution, lbcs)

  soln = mx.Solution()
  soln.loadEigenmode(type, modename, dir)
  soln.grid = mx.Grid(3 * [resolution], origin, [lx, ly, lz])
  return soln.eigFieldValues(type, modename, pointSet, pointSetName, dir)

def analyticFieldValuesAtPoints(dir, modename, pointSet, pointSetName):
  create = False
  # calculate analytic fields at points if data does not already exist in simulation directory
  anFieldFname = "elecField_" + pointSetName + "_" + modename + ".h5"
  if os.path.exists(dir + anFieldFname):
    h5 = tables.openFile(dir + anFieldFname, 'r')
    anField = h5.root.data.read()
    h5.close()
    if len(anField) == len(pointSet):
      return anField

  pol, n, l, m, phase = dm.strToMode(modename)
  anField = dm.efield(pol, n, l, m, phase, pointSet)
  h5 = tables.openFile(dir + anFieldFname, 'w')
  h5.createArray("/", "data", anField)
  h5.close()
  return anField

#def normalizeGrid(resolution, lbcs, modeNames):
  

def normalizeMinRes(resolution, lbcs, modeNames):
  """
  Return a list of coefficients that minimize the L2 error of the 
  bulk electric fields for modeNames
  """

  nModes = len(modeNames)

  volPtsIn = getVolumePoints(32, 64, "in")
  volPtsOut = getVolumePoints(32, 64, "out")
  print "Total number of volume points for minimizing L2 error: ", len(volPtsIn) + len(volPtsOut)

  # get mx fields from simulation
  # calculate analytic fields at points if data does not already exist in simulation directory
  mxVolFieldsIn = []
  mxVolFieldsOut = []
  anVolFieldsIn = []
  anVolFieldsOut = []
  for mode in modeNames:
    print "  Mode " + mode + ":"
    print "  calculating simulation volume field values..."
    mxVolFieldsIn.append(mxFieldValuesAtPoints(resolution, lbcs, 'E', mode, volPtsIn, "volume_in").ravel())
    mxVolFieldsOut.append(mxFieldValuesAtPoints(resolution, lbcs, 'E', mode, volPtsOut, "volume_out").ravel())
    print "  calculating analytic volume field values..."
    anVolFieldsIn.append(analyticFieldValuesAtPoints("./analytic/", mode, volPtsIn, "volume_in").ravel())
    anVolFieldsOut.append(analyticFieldValuesAtPoints("./analytic/", mode, volPtsOut, "volume_out").ravel())
  mxVolFieldsIn = numpy.asarray(mxVolFieldsIn)
  mxVolFieldsOut = numpy.asarray(mxVolFieldsOut)
  anVolFieldsIn = numpy.asarray(anVolFieldsIn)
  anVolFieldsOut = numpy.asarray(anVolFieldsOut)

  # dot product matrix of degenerate mx fields
  dotProds = eps * numpy.dot(mxVolFieldsIn, mxVolFieldsIn.T)
  dotProds += numpy.dot(mxVolFieldsOut, mxVolFieldsOut.T)
  invDotProds = numpy.linalg.inv(dotProds)

  rhs = eps * numpy.dot(mxVolFieldsIn, anVolFieldsIn.T)
  rhs += numpy.dot(mxVolFieldsOut, anVolFieldsOut.T)
  return numpy.dot(invDotProds, rhs)

def eFieldErrorsAtPoints(resolution, lbcs, modeNames, pointSet, pointSetName):
  """
  Return 1,2,Inf-Norms for each mode in modeNames where all modes listed in modeNames
  are assumed to be degenerate. 
  """

  print "\n\nCalculating errors of fields: "
  print "  ", resolution
  print "  ", lbcs
  print "  ", modeNames
  print "  ", pointSetName

  dir = simDir(resolution, lbcs)

  nModes = len(modeNames)

  # get points for finding discrete field coeffs that minimize the L2 volume error
  volPtsIn = getVolumePoints(32, 64, "in")
  volPtsOut = getVolumePoints(32, 64, "out")
  print "Total number of volume points for minimizing L2 error: ", len(volPtsIn) + len(volPtsOut)
  print "Total number of surface points: ", len(pointSet)

  # get mx fields from simulation
  mxSurfFields = []
  mxVolFieldsIn = []
  mxVolFieldsOut = []
  for mode in modeNames:
    print "  Mode " + mode + ":"
    print "  calculating simulation surface field values..."
    mxSurfFields.append(mxFieldValuesAtPoints(resolution, lbcs, 'E', mode, pointSet, pointSetName).ravel())
    print "  calculating simulation volume field values..."
    mxVolFieldsIn.append(mxFieldValuesAtPoints(resolution, lbcs, 'E', mode, volPtsIn, "volume_in").ravel())
    mxVolFieldsOut.append(mxFieldValuesAtPoints(resolution, lbcs, 'E', mode, volPtsOut, "volume_out").ravel())
  mxSurfFields = numpy.asarray(mxSurfFields)

  # calculate analytic fields at points if data does not already exist in simulation directory
  anSurfFields = []
  anVolFieldsIn = []
  anVolFieldsOut = []
  for mode in modeNames:
    print "  Mode " + mode + ":"
    print "  calculating analytic surface field values..."
    anSurfFields.append(analyticFieldValuesAtPoints(dir, mode, pointSet, pointSetName).ravel())
    print "  calculating analytic volume field values..."
    anVolFieldsIn.append(analyticFieldValuesAtPoints("./analytic/", mode, volPtsIn, "volume_in").ravel())
    anVolFieldsOut.append(analyticFieldValuesAtPoints("./analytic/", mode, volPtsOut, "volume_out").ravel())

  # dot product matrix of degenerate mx fields
  dotProds = numpy.zeros((nModes, nModes))
  for i in range(nModes):
    for j in range(i, nModes):
      dotProds[i, j] = numpy.dot(mxVolFieldsIn[i], mxVolFieldsIn[j]) + numpy.dot(mxVolFieldsOut[i], mxVolFieldsOut[j])
      dotProds[j, i] = dotProds[i, j]
  invDotProds = numpy.linalg.inv(dotProds)

  #print normalizeMinRes(resolution, lbcs, modeNames)

  # calculate errors
  norm1s = []
  norm2s = []
  normInfs = []
  for anVolFieldIn, anVolFieldOut, anSurfField in zip(anVolFieldsIn, anVolFieldsOut, anSurfFields):
    rhs = numpy.zeros(nModes)
    for i in range(nModes):
      rhs[i] = numpy.dot(anVolFieldIn, mxVolFieldsIn[i]) + numpy.dot(anVolFieldOut, mxVolFieldsOut[i])
    coeffs = numpy.dot(invDotProds, rhs)
    mxSurfField = numpy.dot(mxSurfFields.T, coeffs)
    diff = mxSurfField - anSurfField
    norm1s.append(numpy.linalg.norm(diff, 1) / numpy.linalg.norm(anSurfField, 1))
    norm2s.append(numpy.linalg.norm(diff, 2) / numpy.linalg.norm(anSurfField, 2))
    normInfs.append(numpy.linalg.norm(diff, numpy.inf) / numpy.linalg.norm(anSurfField, numpy.inf))
  
  return norm1s, norm2s, normInfs

def analyzeResults(resolution, lbcs, pointSets):
  dir = simDir(resolution, lbcs)

  modeNames, modeFreqs = dm.octModes(lbcs, nev)

  # find all degeneracies
  degens = {}
  for modeName in modeNames:
    #polnl = modeName[:4]
    pol, n, l, m, phase = dm.strToMode(modeName)
    polnl = (pol, n, l)
    if degens.has_key(polnl):
      degens[polnl].append(modeName)
    else:
      degens[polnl] = [modeName]

  for ptSetName, points in pointSets.iteritems():
    f = open(dir + "fieldErrors_" + ptSetName + ".txt", 'w')
    for nl, names in degens.iteritems():
      norm1s, norm2s, normInfs = eFieldErrorsAtPoints(resolution, lbcs, names, points, ptSetName)
      for name, norm1, norm2, normInf in zip(names, norm1s, norm2s, normInfs):
        f.write(name + ", %.8e, %.8e, %.8e\n" % (norm1, norm2, normInf))
    f.close()


def runMx(resolution):
  maxSize = 48**3
  size = resolution**3
  np = size / maxSize + 1
  if (np == 1): 
    os.system("../maxwell/build2/src/maxwell --infile=./" + simName + ".mx")
  else:
    #os.system("mpirun -machinefile ./nodes -np " + str(np) + " ../src/maxwell --infile=./" + simName + ".mx")
    os.system("mpirun -machinefile ./nodes -np 2 ../src/maxwell --infile=./" + simName + ".mx")
  
def run(ns):
  #allLBCs = [enumBCs(["pmc", "pec"], 3)[0]]
  #allLBCs = enumBCs(["pmc", "pec"], 3)
  allLBCs = [["pmc", "pmc", "pec"]]

  for n in ns:
    pointSets = {}
    for pad in shellPads:
      pointSets["shell_%.2ddx" % pad] = getShellPoints(n, 32, 100, pad)
    for lbcs in dm.octLBCs:
    #for lbcs in allLBCs:
      writeSim(n, lbcs, pointSets)
      runMx(n)
      saveResults(n, lbcs)
      #analyzeResults(n, lbcs, pointSets)

if __name__ == "__main__":
  #ns = [10, 12, 16, 20]
  #ns = [16, 32, 64, 128]
  ns = [24]
  #ns = [32, 40, 48, 56, 64, 80, 96, 128]
  run(ns)
