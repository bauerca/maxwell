import math
import numpy
from scipy.special import jn, yn
from scipy.optimize import brentq
#import pylab
import sys
import os

epsIn = 10.0
epsOut = 1.0
rtEpsIn = math.sqrt(epsIn)
rtEpsOut = math.sqrt(epsOut)
a = 0.37
b = 0.49

lightspeed = 299792458.

def factorial(n):
  res = 1
  for i in range(1, n + 1):
    res *= i
  return res

def toCyl(pts):
  pts = numpy.asarray(pts)
  rs = numpy.sqrt(numpy.sum(pts[...,0:2]**2, axis = -1))
  phis = numpy.arctan2(pts[...,1], pts[...,0])
  zs = pts[...,2]
  return numpy.array([rs, phis, zs])

def djn(n, x, d):
  if d == 0:
    return jn(n, x)
  else:
    return 0.5 * (djn(n - 1, x, d - 1) - djn(n + 1, x, d - 1))

def dyn(n, x, d):
  if d == 0:
    return yn(n, x)
  else:
    return 0.5 * (dyn(n - 1, x, d - 1) - dyn(n + 1, x, d - 1))

def kFuncTE(k, m):
  kaIn = rtEpsIn * k * a
  kaOut = rtEpsOut * k * a
  kbOut = rtEpsOut * k * b

  one = rtEpsIn * jn(m, kaIn)
  one *= djn(m, kaOut, 1) * dyn(m, kbOut, 1) - djn(m, kbOut, 1) * dyn(m, kaOut, 1)

  two = rtEpsOut * djn(m, kaIn, 1)
  two *= jn(m, kaOut) * dyn(m, kbOut, 1) - djn(m, kbOut, 1) * yn(m, kaOut)

  return one - two

def kFuncTM(k, m):
  kaIn = rtEpsIn * k * a
  kaOut = rtEpsOut * k * a
  kbOut = rtEpsOut * k * b

  one = rtEpsIn * djn(m, kaIn, 1)
  one *= jn(m, kaOut) * yn(m, kbOut) - jn(m, kbOut) * yn(m, kaOut)

  two = rtEpsOut * jn(m, kaIn)
  two *= djn(m, kaOut, 1) * yn(m, kbOut) - jn(m, kbOut) * dyn(m, kaOut, 1)

  return one - two

def findK(pol, m, n):
  #print "Finding " + modeToStr(pol, m, n, "re") + " k-value..."
  if pol == "TE":
    func = kFuncTE
  elif pol == "TM":
    func = kFuncTM

  # asymptotic form for bessel functions goes like a cosine
  kres = 40
  krtgap = math.pi / (b * max([rtEpsIn, rtEpsOut]))
  kstep = krtgap / float(kres)
  kmin = kstep
  kmax = krtgap
  root = 0
  while True:
    #print "searching range: ", kmin, " to ", kmax
    ks = numpy.arange(kmin, kmax, kstep)
    fs = func(ks, m)

    signs = numpy.unique(numpy.sign(fs))
    #print "  signs: ", signs
    roots = int(numpy.sum(numpy.abs(signs))) - 1
    #print "  last root index: ", root
    #print "  range has " + str(roots) + " roots"

    # root is in current range. find it.
    if root + roots >= n:
      kroot = 0.0
      for k1, k2, f1, f2 in zip(ks[:-1], ks[1:], fs[:-1], fs[1:]):
        if f1 * f2 > 0.0:
          continue
        elif f1 == 0.0:
          kroot = k1
        elif f2 == 0.0:
          kroot = k2
        else:
          kroot = brentq(func, k1, k2, args=(m, ))
        root += 1

        if root == n:
          return kroot
    else:
      root += roots

    kmin = kmax
    kmax += krtgap
      

def findKs(pol, kmax, kres = 100, save = True):
  """
  Return all k-values (where k = \omega / c, and \omega is an eigenfrequency) up to 'kmax'
  for a given polarization. The return value is a 2d numpy array of k-values accessed by
  [n, l] where n is the radial 'quantum' number and l is the polar 'quantum' number.
  """

  # approx every pi, there is a root for the spherical bessel functions, where the 
  # argument is sqrt(epsMax) * k * b at the largest
  kRtGap = math.pi / (b * max([rtEpsIn, rtEpsOut]))
  nmax = int(kmax / kRtGap)
  nmax *= 2 # safety factor of 2

  kstep = math.pi / (kRtGap * float(kres))
  print "Taking " + str(kmax / kstep) + " steps per l value"

  # how does the lowest k scale with l?
  # Jackson says for cyl bessel functions: x_vn = n*pi + (v - 1/2)*pi/2
  # let n = 1, and v -> l + 1/2
  # then, x_v1 = (l/2 + 1)*pi
  # sqrt(epsMax)*kmax*b = (lmax/2 + 1)*pi
  # lmax = 2*(sqrt(epsMax)*kmax*b/pi - 1)
  lmax = 2 * int(max([rtEpsIn, rtEpsOut]) * kmax * b / math.pi - 1)
  lmax *= 2 # safety factor of 2

  #ks = []
  ks = numpy.zeros((nmax, lmax), dtype = 'd')

  if pol == "TE":
    func = kFuncTE
  elif pol == "TM":
    func = kFuncTM
  else:
    print "findKs(...): Wrong polarization. Should be 'TE' or 'TM'"

  l = 1
  nextL = True
  while nextL:
    print "l = " + str(l)
    lKs = []
    k1 = kstep
    k2 = k1 + kstep
    step = 0
    n = 0
    while k1 < kmax:
      f1 = func(k1, l)
      f2 = func(k2, l)
      if f1 * f2 < 0.0:
        print "  Step " + str(step) + ": Finding root between k = " + str(k1) + " and k = " + str(k2)
        ks[n, l] = brentq(func, k1, k2, args=(l, ))
        n += 1
        #lKs.append(scipy.optimize.brentq(func, k1, k2, args=(l, )))
      elif f1 == 0.0:
        ks[n, l] = k1
        n += 1
        #lKs.append(k1)
      elif f2 == 0.0:
        ks[n, l] = k2
        n += 1
        #lKs.append(k2)
      k1 = k2
      k2 += kstep
      step += 1
    #ks.append(lKs)

    if n == 0:
      nextL = False

    l += 1

  #rows = max(map(len, ks)) # number of n qnumbers
  #cols = len(ks) + 1 # number l qnumbers
  #res = numpy.zeros((rows, cols))

  #ks = numpy.asarray(ks).T
  if save:
    numpy.savetxt("dsph_msph_ks_a%f_b%f_epsIn%f_epsOut%f.txt" % (a, b, epsIn, epsOut), ks)
  return ks


def aTM(m, k):
  kaIn = rtEpsIn * k * a
  kaOut = rtEpsOut * k * a
  kbOut = rtEpsOut * k * b

  res = jn(m, kaIn)
  res /= jn(m, kaOut) - jn(m, kbOut) * yn(m, kaOut) / yn(m, kbOut)
  return  res

def bTM(m, k):
  kbOut = rtEpsOut * k * b
  return -aTM(m, k) * jn(m, kbOut) / yn(m, kbOut)

def aTE(m, k):
  kaIn = rtEpsIn * k * a
  kaOut = rtEpsOut * k * a
  kbOut = rtEpsOut * k * b

  res = jn(m, kaIn)
  res /= jn(m, kaOut) - djn(m, kbOut, 1) * yn(m, kaOut) / dyn(m, kbOut, 1)
  return res
  
def bTE(m, k):
  kbOut = rtEpsOut * k * b
  return -aTE(m, k) * djn(m, kbOut, 1) / dyn(m, kbOut, 1)

def efieldTM(m, k, pts, phase = "re"):
  pts = numpy.asarray(pts, dtype = 'd')
  r, phi, z = toCyl(pts)
  if phase == "re":
    trig = numpy.cos(m * phi)
  else:
    trig = numpy.sin(m * phi)

  res = numpy.zeros(r.shape + (3, ), dtype = 'd')

  indsIn = r < a
  res[indsIn, 2] = jn(m, rtEpsIn * k * r[indsIn]) * trig[indsIn]

  indsOut = numpy.logical_and(r >= a, r < b)
  xOut = rtEpsOut * k * r[indsOut]
  res[indsOut, 2] = aTM(m, k) * jn(m, xOut) + bTM(m, k) * yn(m, xOut)
  res[indsOut, 2] *= trig[indsOut]

  return res

def efieldTE(m, k, pts, phase = "re"):
  pts = numpy.asarray(pts, dtype = 'd')
  r, phi, z = toCyl(pts)
  if phase == "re":
    rTrig = -numpy.sin(m * phi)
    phiTrig = numpy.cos(m * phi)
  else:
    rTrig = numpy.cos(m * phi)
    phiTrig = numpy.sin(m * phi)

  res = numpy.zeros(r.shape + (3, ), dtype = 'd')

  indsIn = numpy.logical_and(r > 0.0, r < a)
  xIn = rtEpsIn * k * r[indsIn]
  er = m * jn(m, xIn) * rTrig[indsIn] / (epsIn * r[indsIn])
  print numpy.sum(numpy.abs(er))
  ephi = -k * djn(m, xIn, 1) * phiTrig[indsIn] / rtEpsIn
  print numpy.sum(numpy.abs(ephi))
  phiIn = phi[indsIn]
  res[indsIn, 0] = er * numpy.cos(phiIn) - ephi * numpy.sin(phiIn)
  res[indsIn, 1] = er * numpy.sin(phiIn) + ephi * numpy.cos(phiIn)

  indsOut = numpy.logical_and(r >= a, r < b)
  xOut = rtEpsOut * k * r[indsOut]
  er = m * (aTE(m, k) * jn(m, xOut) + bTE(m, k) * yn(m, xOut)) * rTrig[indsOut] / (epsOut * r[indsOut])
  print numpy.sum(numpy.abs(er))
  ephi = -k * (aTE(m, k) * djn(m, xOut, 1) + bTE(m, k) * dyn(m, xOut, 1)) * phiTrig[indsOut] / rtEpsOut
  print numpy.sum(numpy.abs(ephi))
  phiOut = phi[indsOut]
  res[indsOut, 0] = er * numpy.cos(phiOut) - ephi * numpy.sin(phiOut)
  res[indsOut, 1] = er * numpy.sin(phiOut) + ephi * numpy.cos(phiOut)

  return res


def efield(pol, m, n, phase, pts):
  k = findK(pol, m, n)
  if pol == "TE":
    return efieldTE(m, k, pts, phase)
  elif pol == "TM":
    return efieldTM(m, k, pts, phase)
  else:
    print "dsph_in_msph.efield(...): not a valid polarization."


testPt = [0.1, 0.2, 0.3]
#print numpy.cross(testPt, Xlm(2, 1, testPt) + (-1.0)**1 * Xlm(2, -1, testPt))

#print efieldTM(2, 1, 4.3165735561686065791, [0.1, 0.2, 0.3], "re")
#print bfieldTM(2, 1, 4.3165735561686065791, [0.1, 0.2, 0.3])
#print efieldTE(2, 2, 4.0593020005144315490, [0.1, 0.2, 0.3], "im")
#print efieldTE(2, 2, 4.0593020005144315490, [0.1, 0.2, 0.3], "re")

def _quadLBCsTEre(m):
  if m % 2 == 0:
    return ["pec", "pec"]
  else:
    return ["pmc", "pec"]

def _quadLBCsTEim(m):
  if m % 2 == 0:
    return ["pmc", "pmc"]
  else:
    return ["pec", "pmc"]

def quadLBCs(pol, m, n, phase):
  if pol == "TE":
    if phase == "re":
      return _quadLBCsTEre(m)
    else:
      return _quadLBCsTEim(m)
  else:
    if phase == "re":
      return _quadLBCsTEim(m)
    else:
      return _quadLBCsTEre(m)

def quadModeNames(pol, lbcs, numModes, kmax = 10.0):
  all = modeNames(numModes = 2 * numModes)
  res = []

  for modeName in all:
    modepol, m, n, phase = strToMode(modeName)
    if pol == modepol and lbcs == quadLBCs(pol, m, n, phase):
      res.append(modeName)
      if len(res) == numModes:
        break

  return res

def quadModeFreqs(pol, lbcs, numModes):
  allNames, allFreqs = modeFreqs(2 * numModes, 2 * numModes)
  resNames = []
  resFreqs = []

  for modeName, modeFreq in zip(allNames, allFreqs):
    modepol, m, n, phase = strToMode(modeName)
    if pol == modepol and lbcs == quadLBCs(pol, m, n, phase):
      resNames.append(modeName)
      resFreqs.append(modeFreq)
      if len(resNames) == numModes:
        break

  return resNames, resFreqs

def modeFreqs(mmax, nmax):
  ks = []
  names = []
  for m in range(mmax):
    for n in range(1, nmax + 1):
      kTE = findK("TE", m, n)
      kTM = findK("TM", m, n)
      ks += [kTE]
      ks += [kTM]
      names += [modeToStr("TE", m, n, "re")]
      names += [modeToStr("TM", m, n, "re")]
      if m != 0:
        ks += [kTE]
        ks += [kTM]
        names += [modeToStr("TE", m, n, "im")]
        names += [modeToStr("TM", m, n, "im")]
  ks = numpy.asarray(ks)
  names = numpy.asarray(names)
  
  inds = numpy.argsort(ks)
  ks = numpy.take(ks, inds) * 0.5 * lightspeed / math.pi
  names = numpy.take(names, inds)

  return names.tolist(), ks.tolist()

def modeNames(numModes = 10, kmax = None):
  ks = []
  names = []
  for m in range(numModes):
    for n in range(1, numModes + 1):
      kTE = findK("TE", m, n)
      kTM = findK("TM", m, n)
      ks += [kTE]
      ks += [kTM]
      names += [modeToStr("TE", m, n, "re")]
      names += [modeToStr("TM", m, n, "re")]
      if m != 0:
        ks += [kTE]
        ks += [kTM]
        names += [modeToStr("TE", m, n, "im")]
        names += [modeToStr("TM", m, n, "im")]
  ks = numpy.asarray(ks)
  names = numpy.asarray(names)
  
  inds = numpy.argsort(ks)
  ks = numpy.take(ks, inds)
  names = numpy.take(names, inds)

  return names[:numModes].tolist()


def strToMode(str):
  pol = str[:2]
  m = int(str[2])
  n = int(str[3])
  phase = str[-2:]
  return pol, m, n, phase

def modeToStr(pol, m, n, phase):
  return pol + str(m) + str(n) + phase

def modeToVtk(resolution, modename, pathToH5Utils = ""):
  import tables

  delta = 0.5 / float(resolution)

  ex = numpy.zeros(3 * (resolution, ), dtype = 'd')
  ey = numpy.zeros(3 * (resolution, ), dtype = 'd')
  ez = numpy.zeros(3 * (resolution, ), dtype = 'd')

  pol, n, l, m, phase = strToMode(modename)
  print "Saving %s%d%d%d%s field" % (pol, n, l, m, phase)
  #print "Saving TM%d%d%d field" % (n, l, m)
  fx = "elecField_" + modename + "_x.h5"
  fy = "elecField_" + modename + "_y.h5"
  fz = "elecField_" + modename + "_z.h5"
  fvtk = "elecField_" + modename + ".vtk"
  h5x = tables.openFile(fx, 'w')
  h5y = tables.openFile(fy, 'w')
  h5z = tables.openFile(fz, 'w')
  x = 0.0
  for ix in range(resolution):
    y = 0.0
    for iy in range(resolution):
      z = 0.0
      for iz in range(resolution):
        # do stuff
        ex[ix, iy, iz], ey[ix, iy, iz], ez[ix, iy, iz] = efield(pol, n, l, m, phase, [x, y, z])
        z += delta
      y += delta
    x += delta
  h5x.createArray("/", "data", ex)
  h5y.createArray("/", "data", ey)
  h5z.createArray("/", "data", ez)
  h5x.close()
  h5y.close()
  h5z.close()
  # now use h5tovtk to get vtk file
  os.system(pathToH5Utils + "h5tovtk -o %s %s %s %s" % (fvtk, fx, fy, fz))
  os.system("rm " + fx)
  os.system("rm " + fy)
  os.system("rm " + fz)


if __name__ == "__main__":

  print "starting calculation"
  print findKs("TM", 10.0, kres = 100)
  sys.exit()


  fNames = []
  fVals = []
  print "Frequencies:"
  for kind in ["TM", "TE"]:
    for l in range(1, numL):
      for n in range(numN):
        s = "%s%d%dm" % (kind, n, l)
        fNames.append(s)
        if kind == "TE":
          fVals.append(fTEs[n, l])
        elif kind == "TM":
          fVals.append(fTMs[n, l])
        print "  " + s + " = %.8e" % fVals[-1]
  
  inds = numpy.argsort(fVals)
  fVals = numpy.take(fVals, inds)
  fNames = numpy.take(fNames, inds)
  for fVal, fName in zip(fVals, fNames):
    print fName + "  = %.8e" % fVal
  
  print
  #print fTMs, "\n"
  
  print "TE frequencies:"
  for l in range(1, numL):
    for n in range(numN):
      print "  fTE%d%dm = %.8e" % (n, l, fTEs[n, l])
  print
  #print fTEs, "\n"
  
  nx, ny, nz = n = 3 * [16]
  lx, ly, lz = l = 3 * [0.5]
  dx, dy, dz = d = numpy.asarray(l) / numpy.asarray(n, dtype='d')
  ex = numpy.zeros((nx, ny, nz), dtype = 'd')
  ey = numpy.zeros((nx, ny, nz), dtype = 'd')
  ez = numpy.zeros((nx, ny, nz), dtype = 'd')
  
  import tables

  for pol in ["TE", "TM"]:
    for l in range(1, numL):
      for n in range(numN):
        for m in range(l + 1):
          for phase in ["re", "im"]:
            if m == 0 and phase == "im":
              continue
            print "Saving %s%d%d%d%s field" % (pol, n, l, m, phase)
            #print "Saving TM%d%d%d field" % (n, l, m)
            fx = "elecField_%s%d%d%d%s_comp0.h5" % (pol, n, l, m, phase)
            fy = "elecField_%s%d%d%d%s_comp1.h5" % (pol, n, l, m, phase)
            fz = "elecField_%s%d%d%d%s_comp2.h5" % (pol, n, l, m, phase)
            fvtk = "elecField_%s%d%d%d%s.vtk" % (pol, n, l, m, phase)
            h5x = tables.openFile(fx, 'w')
            h5y = tables.openFile(fy, 'w')
            h5z = tables.openFile(fz, 'w')
            x = 0.0
            for ix in range(nx):
              y = 0.0
              for iy in range(ny):
                z = 0.0
                for iz in range(nz):
                  # do stuff
                  ex[ix, iy, iz], ey[ix, iy, iz], ez[ix, iy, iz] = efield(pol, n, l, m, phase, [x, y, z])
                  z += dz
                y += dy
              x += dx
            h5x.createArray("/", "data", ex)
            h5y.createArray("/", "data", ey)
            h5z.createArray("/", "data", ez)
            h5x.flush()
            h5y.flush()
            h5z.flush()
            h5x.close()
            h5y.close()
            h5z.close()
            # now use h5tovtk to get vtk file
            os.system("h5tovtk -o %s %s %s %s" % (fvtk, fx, fy, fz))
  
  
  # get shell of points for 
  #resolution = 128.
  #dx = 0.5 / resolution
  #dTheta = 0.5 * math.pi / 1000.
  #shellRad = a + 4.0 * dx
  #myN = 0
  #myL = 1
  #myM = 0
  #
  #points = []
  #efield = []
  #bfield = []
  #for theta in numpy.arange(dTheta, 0.5 * math.pi, dTheta):
  #  for phi in numpy.arange(0.0, 0.5 * math.pi, dTheta / math.sin(theta)):
  #    pt = toCart([shellRad, theta, phi])
  #    points.append(pt)
  #    print pt
  #    efield.append(efieldTM(myL, myM, kTMs[myL][myN], pt))
  #    bfield.append(bfieldTM(myL, myM, kTMs[myL][myN], pt))
  
  
  
  
  
  
  
  
  #pylab.plot(ks, kFuncTM(2, ks))
  #pylab.plot(ks, f)
  #pylab.plot(ks, sphBesselJ(2, ks, 1))
  #pylab.show()



