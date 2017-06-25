#
#
#  module for creating Maxwell input files
#
#

import os
import sys
import numpy
import math
import copy

havePyTables = False
try:
  import tables
  havePyTables = True
except:
  havePyTables = False


lightspeed = 2.99792458e8
xhat = [1,0,0]
yhat = [0,1,0]
zhat = [0,0,1]
hats = [xhat, yhat, zhat]
zzz = [0,0,0]

def _to3d(lst, pad = 0):
  newlst = 3 * [pad]
  for i in range(len(lst)):
    newlst[i] = lst[i]
  return newlst

def _files(dir, substrs):
  files = os.listdir(dir)
  for s in substrs:
    files = [file for file in files if s in file]
  files.sort()
  return files


def _indent(s, ntabs):
  lines = s.splitlines(True) # keep '\n' at ends
  lines = [ntabs * "  " + line for line in lines]
  return ''.join(lines)


def _xml_string(tag, attrs, children):
  s = "<" + tag + ">\n"
  for name,val in attrs:
    s += "  " + str(name) + " = " + str(val) + "\n"
  for child in children:
    s += _indent(str(child), 1)
  s += "</" + tag + ">\n"
  return s

class Object(object):
  def __init__(self, tag):
    self.tag = str(tag)
    self.attrs = []
    self.children = []

  def _get_attrs(self):
    return []

  def _get_children(self):
    return []

  def __str__(self):
    return _xml_string(self.tag, self._get_attrs(), self._get_children())

class Grid(Object):
  def __init__(self, resolution, origin, dimensions):
    Object.__init__(self, "Grid")
    self.reset(resolution, origin, dimensions)

  def reset(self, resolution, origin, dimensions):
    self.n = list(resolution)
    self.o = list(origin)
    self.l = list(dimensions)
    self.d = map(lambda li, ni: float(li) / float(ni), self.l, self.n)

  def _get_attrs(self):
    attrs = [("resolution", self.n), ("origin", self.o), ("dimensions", self.l)]
    return super(Grid, self)._get_attrs() + attrs


#  def _str(self, dim, ntabs):
#    s = "<Grid>\n"
#    s += "  resolution = " + str(self.n[:dim]) + "\n"
#    s += "  origin = " + str(self.o[:dim]) + "\n"
#    s += "  dimensions = " + str(self.l[:dim]) + "\n"
#    s += "</Grid>\n"
#    return _indent(s, ntabs)

class Dielectric(Object):
  def __init__(self, shape, epsDiag = [1, 1, 1], epsOffDiag = [0, 0, 0], losstan = 0):
    Object.__init__(self, "Dielectric")
    self.epsDiag = epsDiag
    self.epsOffDiag = epsOffDiag
    self.losstan = losstan
    self.shape = shape

  def _get_attrs(self):
    attrs = [("eps diag", list(self.epsDiag)), ("eps off diag", list(self.epsOffDiag)),
      ("loss tangent", float(self.losstan))]
    return super(Dielectric, self)._get_attrs() + attrs

  def _get_children(self):
    children = [self.shape]
    return super(Dielectric, self)._get_children() + children

class Mu(Object):
  def __init__(self, shape, muDiag = [1, 1, 1], muOffDiag = [0, 0, 0], losstan = 0):
    Object.__init__(self, "Mu")
    self.muDiag = muDiag
    self.muOffDiag = muOffDiag
    self.losstan = losstan
    self.shape = shape

  def _get_attrs(self):
    attrs = [("mu diag", list(self.muDiag)), ("mu off diag", list(self.muOffDiag)),
      ("loss tangent", float(self.losstan))]
    return super(Mu, self)._get_attrs() + attrs

  def _get_children(self):
    children = [self.shape]
    return super(Mu, self)._get_children() + children

class PML(Object):
  def __init__(self, shape, alphaMax, exponent, profileBegin, profileEnd):
    super(PML, self).__init__("PML")
    self.shape = shape
    self.alphaMax = alphaMax
    self.exponent = exponent
    self.begin = profileBegin
    self.end = profileEnd

  def _get_attrs(self):
    attrs = [("max alpha", float(self.alphaMax)),
      ("exponent", float(self.exponent)),
      ("profile begin", list(self.begin)),
      ("profile end", list(self.end))]
    return super(PML, self)._get_attrs() + attrs

  def _get_children(self):
    children = [self.shape]
    return super(PML, self)._get_children() + children


def make_pml_box(lowerBounds, upperBounds, profileDir, increasing, alphaMax, exponent):
  """This function returns a simple Cartesian box PML object.

  Parameters:
    lowerBounds : real-space coordinate (list of floats)
    upperBounds : real-space coordinate (list of floats)
    profileDir : integer 0, 1, or 2 for x, y, or z
    increasing : True or False; profile increases toward upper bound
    alphaMax : maximum imaginary part (of relative mu/epsilon)
    exponent : profile exponent

  Returns:
    mx.PML object

  """
  shape = make_box2(lowerBounds, upperBounds)
  dir = int(profileDir)

  begin = copy.deepcopy(zzz)
  begin[dir] = lowerBounds[dir]

  end = copy.deepcopy(zzz)
  end[dir] = upperBounds[dir]

  if increasing:
    return PML(shape, alphaMax, exponent, begin, end)
  else:
    return PML(shape, alphaMax, exponent, end, begin)

def make_boundary_pmls(dim, innerLBs, innerUBs, outerLBs, outerUBs, alphaMax, exponent):
  """Returns a list of PML objects to be passed to the simulation
  object. The resulting objects specify a cubic shell of PMLs."""

  pmls = []

  for i in range(dim):
    norm = copy.deepcopy(zzz)
    norm[i] = 1.0

    # lower
    lb = outerLBs[i]
    ub = innerLBs[i]
    lbvec = copy.deepcopy(zzz)
    lbvec[i] = lb
    ubvec = copy.deepcopy(zzz)
    ubvec[i] = ub
    midPt = copy.deepcopy(zzz)
    midPt[i] = 0.5*(lb + ub)
    shape = Slab(ub - lb, norm, midPt)
    pmls.append(PML(shape, alphaMax, exponent, ubvec, lbvec))

    # upper
    lb = innerUBs[i]
    ub = outerUBs[i]
    lbvec = copy.deepcopy(zzz)
    lbvec[i] = lb
    ubvec = copy.deepcopy(zzz)
    ubvec[i] = ub
    midPt = copy.deepcopy(zzz)
    midPt[i] = 0.5*(lb + ub)
    shape = Slab(ub - lb, norm, midPt)
    pmls.append(PML(shape, alphaMax, exponent, lbvec, ubvec))

  return pmls

# Geometric transformations section
class Rotation(Object):
  def __init__(self, axis, angle, pivot = None):
    super(Rotation, self).__init__("Rotate")
    self.axis = list(axis)
    self.angle = float(angle)
    if pivot is not None:
      self.pivot = list(pivot)

  def _get_attrs(self):
    attrs = [("type", "rotation"), ("axis", self.axis),
      ("angle", self.angle)]
    return super(Rotation, self)._get_attrs() + attrs

class Translation(Object):
  def __init__(self, vector):
    super(Translation, self).__init__("Translate")
    self.vector = vector

  def _get_attrs(self):
    attrs = [("type", "translation"), ("vector", list(self.vector))]
    return super(Translation, self)._get_attrs() + attrs

class Scale(Object):
  def __init__(self, magnitudes):
    super(Scale, self).__init__("Scale")
    self.magnitudes = magnitudes

  def _get_attrs(self):
    attrs = [("type", "scaling"), ("magnitudes", list(self.magnitudes))]
    return super(Scale, self)._get_attrs() + attrs


class Reflection(Object):
  def __init__(self, normal, pointInPlane):
    super(Reflection, self).__init__("Reflect")
    self.normal = normal
    self.pointInPlane = pointInPlane

  def _get_attrs(self):
    attrs = [("type", "reflection"), ("normal", list(self.normal)),
      ("point in plane", list(self.pointInPlane))]
    return super(Reflection, self)._get_attrs() + attrs



class Shape(Object):
  def __init__(self, tag, name):
    super(Shape, self).__init__(tag)

    self.transforms = []
    self.name = str(name)
    self.inv = "false"

  def invert(self):
    if self.inv == "true":
      self.inv = "false"
    else:
      self.inv = "true"

  def rotate(self, axis, angle, pivot = None):
    self.transforms.append(Rotation(axis, angle, pivot))

  def translate(self, vector):
    self.transforms.append(Translation(vector))

  def scale(self, magnitudes):
    self.transforms.append(Scale(magnitudes))

  def reflect(self, normal, pointInPlane):
    self.transforms.append(Reflection(normal, pointInPlane))

  def _get_attrs(self):
    attrs = [("name", self.name), ("invert", self.inv)]
    return super(Shape, self)._get_attrs() + attrs

  def _get_children(self):
    return super(Shape, self)._get_children() + self.transforms

class Sphere(Shape):
  def __init__(self, radius, location, name = "sphere"):
    super(Sphere, self).__init__("Sphere", name)
    self.radius = radius
    self.translate(location)

  def _get_attrs(self):
    return [("radius", float(self.radius))] + super(Sphere, self)._get_attrs()


class Cylinder(Shape):
  def __init__(self, radius, axis, location, name = "cylinder"):
    super(Cylinder, self).__init__("Cylinder", name)
    self.radius = radius
    self.axis = axis
    self.translate(location)

  def _get_attrs(self):
    attrs = [("radius", float(self.radius)), ("axis", list(self.axis))]
    return attrs + super(Cylinder, self)._get_attrs()


class Torus(Shape):
  def __init__(self, majRadius, minRadius, axis, location, name = "torus"):
    super(Torus, self).__init__("Torus", name)
    self.maj_radius = majRadius
    self.min_radius = minRadius
    self.axis = axis
    self.translate(location)

  def _get_attrs(self):
    attrs = [("minor radius", float(self.min_radius)),
             ("major radius", float(self.maj_radius)),
             ("axis", list(self.axis))]
    return attrs + super(Torus, self)._get_attrs()


class Cone(Shape):
  def __init__(self, angle, axis, vertex, name = "cone"):
    super(Cone, self).__init__("Cone", name)
    self.angle = angle
    self.axis = axis
    self.translate(vertex)

  def _get_attrs(self):
    attrs = [("angle", float(self.angle)), ("axis", list(self.axis))]
    return attrs + super(Cone, self)._get_attrs()

class HalfSpace(Shape):
  def __init__(self, normal, pointInPlane, name = "plane"):
    super(HalfSpace, self).__init__("HalfSpace", name)
    self.normal = normal
    self.translate(pointInPlane)

  def _get_attrs(self):
    return [("normal", list(self.normal))] + super(HalfSpace, self)._get_attrs()

class Slab(Shape):
  def __init__(self, thickness, normal, midplanePoint, name = "slab"):
    super(Slab, self).__init__("Slab", name)
    self.thickness = thickness
    self.normal = normal
    self.translate(midplanePoint)

  def _get_attrs(self):
    attrs = [("thickness", float(self.thickness)), ("normal", list(self.normal))]
    return attrs + super(Slab, self)._get_attrs()

def make_box(sides, center, name="box"):
  slabs = []
  for side,xi,hat in zip(sides, center, hats):
    slabs.append(Slab(side, hat, xi))
  shape = ShapeIntersection(slabs, name=name)
  return shape

def make_box2(lowerBounds, upperBounds, name="box"):
  slabs = []
  for lbi, ubi, hat in zip(lowerBounds, upperBounds, hats):
    midPt = 0.5 * (lbi + ubi)
    width = ubi - lbi
    slabs.append(Slab(width, hat, midPt))
  shape = ShapeIntersection(slabs, name=name)
  return shape

class ShapeSubtract(Shape):
  def __init__(self, base, subtract, name = "subtract"):
    super(ShapeSubtract, self).__init__("ShapeSubtract", name)
    self.base = base
    self.subtract = subtract

  def setBase(self, shape):
    self.base = shape

  def setSubtraction(self, shape):
    self.subtract = shape

  def _get_children(self):
    # put transform objects after child shapes
    return [self.base, self.subtract] + super(ShapeSubtract, self)._get_children()

class ShapeUnion(Shape):
  def __init__(self, shapes = [], name = "union"):
    super(ShapeUnion, self).__init__("ShapeUnion", name)
    self.shapes = shapes

  def add(self, shape):
    self.shapes.append(shape);

  def _get_children(self):
    return self.shapes + super(ShapeUnion, self)._get_children()

class ShapeIntersection(Shape):
  def __init__(self, shapes = [], name = "intersection"):
    super(ShapeIntersection, self).__init__("ShapeIntersection", name)
    self.shapes = shapes

  def add(self, shape):
    self.shapes.append(shape);

  def _get_children(self):
    return self.shapes + super(ShapeIntersection, self)._get_children()


class ShapeMirror(Shape):
  def __init__(self, shape, normal, pointInPlane, name = "mirror"):
    super(ShapeMirror, self).__init__("ShapeMirror", name)
    self.n = normal
    self.p = pointInPlane
    self.shape = shape

  def _get_attrs(self):
    attrs = [("normal", list(self.n)),
             ("point in plane", list(self.p))]
    return attrs + super(ShapeMirror, self)._get_attrs()

  def _get_children(self):
    return [self.shape] + super(ShapeMirror, self)._get_children()

class ShapeRepeat(Shape):
  def __init__(self, shape, origin, direction, step, npos, nneg, name = "repeat"):
    super(ShapeRepeat, self).__init__("ShapeRepeat", name)
    self.shape = shape
    self.origin = list(origin)
    self.direction = list(direction)
    self.step = float(step)
    self.npos = int(npos)
    self.nneg = int(nneg)

  def _get_attrs(self):
    attrs = [("origin", list(self.origin)),
             ("direction", list(self.direction)),
             ("step", float(self.step)),
             ("pos steps", int(self.npos)),
             ("neg steps", int(self.nneg))]
    return attrs + super(ShapeRepeat, self)._get_attrs()


  def _get_children(self):
    return [self.shape] + super(ShapeRepeat, self)._get_children()

class PEC(Object):
  def __init__(self, shape, dmFrac=0.0):
    super(PEC, self).__init__("PEC")
    self.shape = shape
    self.dmFrac = dmFrac

  def _get_attrs(self):
    attrs = [("D-M frac", float(self.dmFrac))]
    return attrs + super(PEC, self)._get_attrs()

  def _get_children(self):
    return [self.shape] + super(PEC, self)._get_children()

class BoundaryConditions(Object):
  def __init__(self):
    super(BoundaryConditions, self).__init__("BoundaryConditions")
    self.options = ["periodic", "pec", "pmc"]
    self.uBCs = [ "periodic", "periodic", "periodic"]
    self.lBCs = [ "periodic", "periodic", "periodic"]
    self.phaseShifts = [0.0, 0.0, 0.0]

  def setUpperBCs(self, bcs):
    for bc in bcs:
      if bc not in self.options:
        print "Boundary condition, '" + bc + "', not a valid option.\n"
        sys.exit()
    self.uBCs = bcs

  def setLowerBCs(self, bcs):
    for bc in bcs:
      if bc not in self.options:
        print "Boundary condition, '" + bc + "', not a valid option.\n"
        sys.exit()
    self.lBCs = bcs

  def setPhaseShifts(self, phaseShifts):
    self.phaseShifts = phaseShifts

  def _get_attrs(self):
    attrs = [("lower bcs", "[ " + ", ".join(self.lBCs) + "]"),
             ("upper bcs", "[ " + ", ".join(self.uBCs) + "]"),
             ("phase shifts", list(self.phaseShifts))]
    return attrs + super(BoundaryConditions, self)._get_attrs()


class Output(Object):
  def __init__(self):
    super(Output, self).__init__("Output")
    self.directory = "./"
    self.pointSets = {}
    self.outputFields = "true"
    self.outputFreqs = "true"

  def _get_attrs(self):
    attrs = [("fields", self.outputFields),
             ("frequencies", self.outputFreqs)]
    return attrs + super(Output, self)._get_attrs()

#  def fieldValues(self, name, points):
#    self.pointSets[name] = points
#
#  def _str(self, dim, ntabs):
#    tab = "  "
#    s = "<Output>\n"
#    s += "  fields = %s\n" % self.outputFields
#    s += "  frequencies = %s\n" % self.outputFreqs
#    for name, pointSet in self.pointSets.iteritems():
#      s += "  <FieldValues>\n"
#      s += "    name = %s\n" % name
#      if pointSet.__class__ == Grid:
#        s += "    kind = grid\n"
#        s += pointSet._str(dim, 2)
#      else:
#        s += "    kind = points\n"
#        s += 2 * tab + "<Points>\n"
#        p = 0
#        for point in pointSet:
#          s += 3 * tab + "r%d = %s\n" % (p, str(list(point)))
#          p += 1
#        s += 2 * tab + "</Points>\n"
#      s += "  </FieldValues>\n"
#    s += "</Output>\n"
#    return _indent(s, ntabs)

class Eigensolver(Object):
  def __init__(self, params = None, linearSolver = None):
    super(Eigensolver, self).__init__("Eigensolver")
    if params == None:
      self.params = {"type": "krylov-schur",
                      "invert": "true",
                      "nev": 12,
                      "block size": 1,
                      "max restarts": 100,
                      "basis": 20,
                      "tol": 1.0e-6,
                      "spectrum": "LM",
                      "output": 2}
    else:
      self.params = params

    if linearSolver == None:
      self.linearSolver = LinearSolver()
      self.linearSolver.prec = AMG()
    else:
      self.linearSolver = linearSolver

  def setParams(self, params):
    for key,val in params.iteritems():
      self.params[key] = val

  def setLinearSolver(self, solver):
    self.linearSolver = solver

  def _get_attrs(self):
    attrs = [pair for pair in self.params.iteritems()]
    return attrs + super(Eigensolver, self)._get_attrs()

  def _get_children(self):
    return [self.linearSolver] + super(Eigensolver, self)._get_children()



class LinearSolver(Object):
  def __init__(self, solveType="gmres", basis=40, tol=1.0e-6, prec=None):
    """
    Parameters:
      solveType : options are 'gmres', 'bicgstab'
      basis : size of Krylov basis used in iterative solve
      tol : convergence tolerance
      prec : a Preconditioner object

    """
    super(LinearSolver, self).__init__("LinearSolver")
    self.solveType = solveType
    self.basis = basis
    self.tol = tol
    self.prec = prec

  def setPreconditioner(self, prec):
    self.prec = prec

  def _get_attrs(self):
    attrs = [("type", str(self.solveType)),
             ("basis", int(self.basis)),
             ("tol", float(self.tol))]
    return attrs + super(LinearSolver, self)._get_attrs()

  def _get_children(self):
    c = []
    if self.prec != None:
      c.append(self.prec)

    return c + super(LinearSolver, self)._get_children()



class Preconditioner(Object):
  def __init__(self, precType):
    super(Preconditioner, self).__init__("Preconditioner")
    self.precType = precType

  def _get_attrs(self):
    return [("type", str(self.precType))] + super(Preconditioner, self)._get_attrs()

class AMG(Preconditioner):
  def __init__(self, sweeps=1, smoother="Chebyshev", coarseSmoother="Chebyshev",
      maxLevels=10):
    super(AMG, self).__init__("amg")
    self.sweeps = sweeps
    self.smoother = smoother
    self.coarseSmoother = coarseSmoother
    self.maxLevels = maxLevels

  def _get_attrs(self):
    attrs = [("sweeps", int(self.sweeps)),
             ("smoother", str(self.smoother)),
             ("coarse smoother", str(self.coarseSmoother)),
             ("max levels", int(self.maxLevels))]
    return super(AMG, self)._get_attrs() + attrs

class ILU(Preconditioner):
  def __init__(self, dropTol=1.0e-3, fill=10.0):
    super(ILU, self).__init__("ilu")
    self.dropTol = dropTol
    self.fill = fill

  def _get_attrs(self):
    attrs = [("drop tol", float(self.dropTol)),
             ("fill", float(self.fill))]
    return super(ILU, self)._get_attrs() + attrs

class Simulation(Object):
  def __init__(self, name, dim, pol = "TE"):
    super(Simulation, self).__init__("Maxwell")
    self.name = name
    self.dim = dim
    self.pol = pol
    self.diels = []
    self.mus = []
    self.pmls = []
    self.pec = None
    self.grid = None
    self.bcs = BoundaryConditions()
    self.output = Output()
    self.solver = Eigensolver()

  def setSolver(self, solver):
    self.solver = solver

  def setGrid(self, grid):
    self.grid = grid

  def addDielectric(self, diel):
    self.diels.append(diel)

  def addMu(self, mu):
    self.mus.append(mu)

  def addPMLs(self, pmls):
    """Add a PML object(s) to the simulation.

    Parameters:
      pmls : a list of PML objects or a single PML object

    """
    try:
      self.pmls += list(pmls)
    except:
      self.pmls += [pmls]

  def addPEC(self, pec):
    self.pec = pec

  def setBCs(self, bcs):
    self.bcs = bcs

  def setOutput(self, output):
    self.output = output

  def _get_attrs(self):
    attrs = [("name", self.name), ("dim", int(self.dim)),
             ("polarization", self.pol), ("domain", "frequency")]
    return attrs + super(Simulation, self)._get_attrs()

  def _get_children(self):
    children = []
    if self.grid == None:
      sys.exit("Simulation needs a Grid")
    else:
      children += [self.grid]

    if self.pec:
      children += [self.pec]

    children += self.diels + self.mus + self.pmls

    if not self.solver.params.has_key("shift"):
      self.solver.params["shift"] = 0.05 * (2.0 * math.pi / max(self.grid.l))**2
    children += [self.bcs, self.solver, self.output]

    return children + super(Simulation, self)._get_attrs()

  def write(self):
    f = open(self.name + ".mx", 'w')
    f.write(str(self))
    f.close()



def h5tovtk(dir, pathToH5Utils = ""):
  import os

  files = os.listdir(dir)
  #elecfiles = [file for file in files if ("mxYeeElecField" in file and ".h5" in file)]
  #magfiles = [file for file in files if ("mxYeeMagField" in file and ".h5" in file)]
  elecfiles = _files(files, ["mxYeeElecField", ".h5"])
  magfiles = _files(files, ["mxYeeMagField", ".h5"])

  elecfiles.sort()
  magfiles.sort()

  modes = len(elecfiles)

  for i in range(modes):
    print "Mode %.2d:" % i
    for files in [elecfiles, magfiles]:
      file = dir + files[i]
      vtkName = file.rstrip(".h5") + ".vtk"
      #evec = dir + "mxYeeElecField_eigenvectors_vec%.2d.vtk" % i
      print "  creating " + vtkName
      cmd = pathToH5Utils + "h5tovtk -o " + vtkName
      cmd += " -t 0 " + file
      cmd += " -t 1 " + file
      cmd += " -t 2 " + file
      os.system(cmd)


class Solution:
  def __init__(self, dir = None):
    self.grid = None
    if dir:
      load(self, dir)

    self._fieldPrefixes = {}
    self._fieldPrefixes['E'] = "mxYeeFitEField"
    self._fieldPrefixes['B'] = "mxYeeFitBField"
    self._fieldPrefixes['D'] = "mxYeeFitDField"
    self._fieldPrefixes['H'] = "mxYeeFitHField"

    self._compLocs = {}
    self._compLocs['E'] = numpy.array([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]])
    self._compLocs['B'] = numpy.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    self._compLocs['D'] = numpy.array([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]])
    self._compLocs['H'] = numpy.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])

    self.eigmodes = {}
    self.eigmodes['E'] = {}
    self.eigmodes['B'] = {}
    self.eigmodes['D'] = {}
    self.eigmodes['H'] = {}

  def loadEigenmode(self, type, modename, dir):
    fieldPrefix = self._fieldPrefixes[type]
    if not self.eigmodes[type].has_key(modename):
      h5 = tables.openFile(dir + "/" + fieldPrefix + "_eigVecs_" + modename + ".h5", "r")
      self.eigmodes[type][modename] = h5.root.real.read()
      h5.close()
    else:
      pass

  def _eigFieldValue(self, type, modename, pt, comp):
    pt = numpy.asarray(_to3d(pt))

    eigmode = self.eigmodes[type][modename]

    o = numpy.asarray(_to3d(self.grid.o))
    d = numpy.asarray(_to3d(self.grid.d))
    n = numpy.asarray(_to3d(self.grid.n, pad = 1))
    l = numpy.asarray(_to3d(self.grid.l))
    cell = (pt - (o + d * self._compLocs[type][comp])) / d
    cell = numpy.floor(cell)
    if numpy.any(cell[:eigmode.ndim - 1] < 0):
      print "point: ", pt
      print "cell: ", cell
      print "gridcell: ", d
      print "origin: ", o
      #exit("Interpolation point is out of simulation bounds")
      print "Interpolation point is out of simulation bounds, returning 0"
      return 0.

    shift = (pt - (cell + self._compLocs[type][comp]) * d) / d
    t, u, v = shift.tolist()

    res = 0.

    if n[2] == 1:
      if eigmode.shape[-1] == 1:
        if comp != 2:
          return 0.0
        else:
          comp = 0
      elif eigmode.shape[-1] == 2:
        if comp == 2:
          return 0.0
      res += (1. - t) * (1. - u) * eigmode[cell[0], cell[1], comp]
      res += t * (1. - u) * eigmode[cell[0] + 1, cell[1], comp]
      res += (1. - t) * u * eigmode[cell[0], cell[1] + 1, comp]
      res += t * u * eigmode[cell[0] + 1, cell[1] + 1, comp]
      return res;
    else:
      res += (1. - t) * (1. - u) * (1. - v) * eigmode[cell[0], cell[1], cell[2], comp]
      res += t * (1. - u) * (1. - v) * eigmode[cell[0] + 1, cell[1], cell[2], comp]
      res += (1. - t) * u * (1. - v) * eigmode[cell[0], cell[1] + 1, cell[2], comp]
      res += (1. - t) * (1. - u) * v * eigmode[cell[0], cell[1], cell[2] + 1, comp]
      res += (1. - t) * u * v * eigmode[cell[0], cell[1] + 1, cell[2] + 1, comp]
      res += t * (1. - u) * v * eigmode[cell[0] + 1, cell[1], cell[2] + 1, comp]
      res += t * u * (1. - v) * eigmode[cell[0] + 1, cell[1] + 1, cell[2], comp]
      res += t * u * v * eigmode[cell[0] + 1, cell[1] + 1, cell[2] + 1, comp]
      return res


  def eigFieldValue(self, type, modename, pt):
    if self.grid == None:
      print "Specify a grid first"
      raise TypeError

    dataDims = list(numpy.shape(self.eigmodes[type][modename]))
    nplus = map(lambda x: x + 1, self.grid.n)[:len(dataDims) - 1]
    if nplus != dataDims[:-1]:
      print "Grid does not match dimensions of data"
      print "  grid n: ", nplus
      print "  data dims: ", dataDims
      raise ValueError

    return map(lambda x: self._eigFieldValue(type, modename, pt, x), range(3))

  def eigFieldValues(self, type, modename, pts, ptsName = None, saveDir = None):
    if saveDir != None:
      if saveDir[-1] != "/":
        saveDir += "/"

    res = []
    # if user supplies name and directory, search for existing data first. Existing
    # data must have the same dimensions of 'pts'
    if saveDir != None and ptsName != None and havePyTables:
      filepath = saveDir + self._fieldPrefixes[type] + "_" + ptsName + "_eigVecs_" + modename + ".h5"
      if os.path.exists(filepath):
        h5 = tables.openFile(filepath, 'r')
        data = h5.root.data.read()
        h5.close()
        if len(data) == len(pts):
          return data

    for pt in pts:
      res.append(self.eigFieldValue(type, modename, pt))
    res = numpy.asarray(res)

    # save if name and directory given
    if saveDir != None and ptsName != None and havePyTables:
      h5 = tables.openFile(filepath, 'w')
      h5.createArray("/", "data", res)
      h5.close()

    return res

  def eigToVtk(self, type, modename, pathToH5Utils = "", dir = "./"):
    files = []
    for comp in range(numpy.shape(self.eigmodes[type][modename])[-1]):
      file = self._fieldPrefixes[type] + "_eigenvectors_" + modename + "_" + str(comp) + ".h5"
      files.append(file)
      h5 = tables.openFile(file, 'w')
      h5.createArray("/", "data", self.eigmodes[type][modename].T[comp].T)
      h5.close()

    vtkName = dir + self._fieldPrefixes[type] + "_eigenvectors_" + modename + ".vtk"
    print "Creating " + vtkName
    os.system(pathToH5Utils + "h5tovtk -o " + vtkName + " " + " ".join(files))

    for file in files:
      os.system("rm " + file)


if __name__ == "__main__":
  sh = ShapeSubtract(Sphere(3, [1, 1, 1]), Sphere(3, [1.2, 1, 1]))
  diel = Dielectric(shape = sh)
  sim = Simulation("phc.mx", 3)
  sim.setGrid(Grid([16, 16, 16], [0, 0, 0], [1, 1, 1]))
  sim.addDielectric(diel)
  sim.write()

  print diel._str(2)
