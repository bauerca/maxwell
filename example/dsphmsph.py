import math
import numpy
from scipy.special import sph_jn, sph_yn, sph_jnyn, lpmv, sph_harm, jv, yv
from scipy.optimize import brentq
#import pylab
import sys
import os

lightspeed = 299792458.

epsIn = 10.0
epsOut = 1.0
rtEpsIn = math.sqrt(epsIn)
rtEpsOut = math.sqrt(epsOut)
a = 0.37
b = 0.49

def factorial(n):
  res = 1
  for i in range(1, n + 1):
    res *= i
  return res

def toSph(pts):
  pts = numpy.asarray(pts)
  rs = numpy.sqrt(numpy.sum(pts**2, axis = -1))
  thetas = numpy.arccos(pts[...,2] / rs)
  phis = numpy.arctan2(pts[...,1], pts[...,0])
  return numpy.array([rs, thetas, phis])

def jl(l, x):
  return numpy.sqrt(0.5 * numpy.pi / x) * jv(l + 0.5, x)

def djl(l, x):
  return jl(l - 1, x) - float(l + 1) * jl(l, x) / x

def yl(l, x):
  return numpy.sqrt(0.5 * numpy.pi / x) * yv(l + 0.5, x)

def dyl(l, x):
  return yl(l - 1, x) - float(l + 1) * yl(l, x) / x

def alt_jl(l, x):
  return x * djl(l, x) + jl(l, x)

def alt_yl(l, x):
  return x * dyl(l, x) + yl(l, x)

def alt_sph_jn(l, x):
  jn, djn = sph_jn(l, x)
  return x * djn[-1] + jn[-1]

def alt_sph_yn(l, x):
  yn, dyn = sph_yn(l, x)
  return x * dyn[-1] + yn[-1]

def kFuncTE(k, l):
  kaIn = rtEpsIn * k * a
  kaOut = rtEpsOut * k * a
  kbOut = rtEpsOut * k * b

  one = rtEpsIn * djl(l, kaIn)
  one = one * (jl(l, kaOut) * yl(l, kbOut) - jl(l, kbOut) * yl(l, kaOut))

  two = rtEpsOut * jl(l, kaIn)
  two = two * (djl(l, kaOut) * yl(l, kbOut) - jl(l, kbOut) * dyl(l, kaOut))

  return one - two

def kFuncTM(k, l):
  kaIn = rtEpsIn * k * a
  kaOut = rtEpsOut * k * a
  kbOut = rtEpsOut * k * b

  one = epsIn * jl(l, kaIn) 
  one *= alt_jl(l, kaOut) * alt_yl(l, kbOut) - alt_jl(l, kbOut) * alt_yl(l, kaOut)

  two = alt_jl(l, kaIn)
  two *= jl(l, kaOut) * alt_yl(l, kbOut) - alt_jl(l, kbOut) * yl(l, kaOut)

  return one - two



def findK(pol, n, l, kstep = 0.01):
  #print "Finding " + modeToStr(pol, m, n, "re") + " k-value..."
  if pol == "TE":
    func = kFuncTE
  elif pol == "TM":
    func = kFuncTM

  krange = math.pi / max([rtEpsIn, rtEpsOut]) / b
  kmin = kstep
  kmax = krange
  root = 0 
  while True:
    #print "searching range: ", kmin, " to ", kmax
    ks = numpy.arange(kmin, kmax, kstep)

    linds, uinds = _rootInds(lambda x: func(x, l), ks)
    roots = len(linds)

    # if root is in current range, find it.
    if root + roots > n:
      lind = linds[n - root]
      uind = uinds[n - root]
      return brentq(func, ks[lind], ks[uind], args = (l, ))
    else:
      root += roots
      kmin = kmax
      kmax += krange


def _rootInds(func, xs):
  """given a range of x-values, return the indices of the x-values in the array
  that bracket a root
  
  returns pair of lists: ([lower indices], [upper indices])
  """
  fs = func(xs)
  signs = numpy.sign(fs)

  if numpy.any(signs == 0.0):
    print "Hit a root, have to code for this now."

  # whenever the signs change, we have a root
  changes = signs[:-1] * signs[1:]
  linds = numpy.nonzero(changes == -1.0)[0]
  uinds = linds + 1

  return linds, uinds


def _kPolInRange(pol, kmin, kmax, kstep = 0.01):
  """
  returns two lists. First is list of modenames, next is list of k-values, both sorted by
  increasing k, of course. This function does not take into acct degeneracies. Format for 
  returned modenames is either 'TM-n-l' or 'TE-n-l'.
  """

  if pol == "TE":
    func = kFuncTE
  elif pol == "TM":
    func = kFuncTM

  allKs = []
  allNames = []

  ksPrev = numpy.arange(kstep, kmin, kstep)
  if kmin == 0.0:
    ks = numpy.arange(kstep, kmax, kstep)
  else:
    ks = numpy.arange(kmin, kmax, kstep)

  l = 1
  while True:
    lKs = []
    lNames = []

    # first find out how many roots have been skipped over for this l number
    n = len(_rootInds(lambda x: func(x, l), ksPrev)[0])

    # next, get bracketing indices for roots in requested range
    linds, uinds = _rootInds(lambda x: func(x, l), ks)

    # we quit looping over l when there are no roots below kmax
    if n + len(linds) == 0:
      inds = numpy.argsort(allKs)
      allKs = numpy.take(allKs, inds)
      allNames = numpy.take(allNames, inds)
      return allNames, allKs
    else:
      for lind, uind in zip(linds, uinds):
        lKs += [brentq(func, ks[lind], ks[uind], args = (l, ))]
        lNames += [pol + "-" + str(n) + "-" + str(l)]
        n += 1

      #print "k-values for l = " + str(l), lKs
      allKs += lKs
      allNames += lNames

      l += 1


def _kInRange(kmin, kmax, kstep = 0.01):
  """
  Same as above except combines both polarizations.
  """

  nameTM, kTM = _kPolInRange("TM", kmin, kmax, kstep)
  nameTE, kTE = _kPolInRange("TE", kmin, kmax, kstep)
  
  k = numpy.concatenate([kTM, kTE])
  names = numpy.concatenate([nameTM, nameTE])

  inds = numpy.argsort(k)
  k = numpy.take(k, inds)
  names = numpy.take(names, inds)

  return names, k


octLBCs = [["pec", "pec", "pec"],
	   ["pec", "pec", "pmc"],
	   ["pec", "pmc", "pec"],
	   ["pmc", "pec", "pec"],
	   ["pec", "pmc", "pmc"],
	   ["pmc", "pec", "pmc"],
	   ["pmc", "pmc", "pec"],
	   ["pmc", "pmc", "pmc"]]

def octModes(lbcs, num):
  """
  Given the -x,-y,-z simulation boundary conditions, this function returns the names and
  frequencies of the lowest 'num' modes that satisfy the bcs, sorted by increasing frequency.
  This function will return more than 'num' modes if returning exactly 'num' modes means
  splitting a group of degenerate modes. 'num' is a lower bound. 
  
  lbcs is a list of 3 boundary conditions, where each b.c. is either 'pmc' or 'pec'
  (example: ['pec', 'pec', 'pmc']). 'pec' is perfect electrical conductor and 'pmc' is
  perfect magnetic conductor.
  """
  
  # first, we have to guess at a kmax that will return 'num' modes that satisfy
  # the boundary conditions

  # the limiting form of the spherical bessel functions is a cosine, so
  # int(rtEpsMax kmax rmax / math.pi) could be a VERY rough approx to the number
  # of modes in the range (only for lowest l number, though)

  krange = float(num) * math.pi / max([rtEpsIn, rtEpsOut]) / b

  names = []
  ks = []

  k1 = 0.0
  k2 = krange
  while len(names) < num:
    tnames, tks = _kInRange(k1, k2)
    for tname, tk in zip(tnames, tks):
      pol = tname.split("-")[0]
      n = int(tname.split("-")[1])
      l = int(tname.split("-")[2])
      for m in range(l + 1):
        for phase in ["re", "im"]:
          if m == 0 and phase == "im":
            continue
          if lbcs == octLowerBCs(pol, l, m, phase):
            names += [tname + "-" + str(m) + "-" + phase]
            ks += [tk]
      if len(names) >= num:
        return names, 0.5 * numpy.array(ks) * lightspeed / math.pi
    k1 = k2
    k2 += krange

  return names, 0.5 * numpy.array(ks) * lightspeed / math.pi


def gradYlm(l, m, pts):
  pts = numpy.asarray(pts)
  r, theta, phi = toSph(pts)

  sinTh = numpy.sin(theta)
  cosTh = numpy.cos(theta)
  sinPh = numpy.sin(phi)
  cosPh = numpy.cos(phi)

  Nlm = math.sqrt(float((2 * l + 1) * factorial(l - m)) / (4.0 * math.pi * float(factorial(l + m))))

  #print "m = ", m, " l = ", l
  # ============
  # The following lines are just for calculating Plm and dPlm, since stupid
  # scipy won't take m < 0, even though -l <= m <= l are valid
  mp = abs(m)
  factor = 1.0
  if m < 0:
    factor = float((-1)**mp * factorial(l - mp)) / float(factorial(l + mp))

  plm = lpmv(mp, l, cosTh)

  if mp == l:
    dplm = -(l + mp) * (l - mp + 1) * numpy.sqrt(1.0 - cosTh**2) * lpmv(mp - 1, l, cosTh)
    dplm -= mp * cosTh * plm
  else:
    dplm = l * cosTh * plm - (l + mp) * lpmv(mp, l - 1, cosTh)
  dplm *= factor / (cosTh**2 - 1.0)
  plm *= factor
  # done calculating Plm and dPlm
  # ==============

  res = numpy.zeros(r.shape + (3, ), dtype = 'complex')

  res[...,0] = -sinTh * cosTh * cosPh * dplm / r - (1.0j) * m * sinPh * plm / (r * sinTh)
  res[...,1] = -sinTh * cosTh * sinPh * dplm / r + (1.0j) * m * cosPh * plm / (r * sinTh)
  res[...,2] = sinTh * sinTh * dplm / r

  return Nlm * numpy.exp(1.0j * m * phi)[...,numpy.newaxis] * res

#print "gradYlm(l=2, m=2, r=[0.1, 0.2, 0.3]): ", gradYlm(2, 2, [0.1, 0.2, 0.3])
#print "gradYlm(l=2, m=-2, r=[0.1, 0.2, 0.3]): ", gradYlm(2, -2, [0.1, 0.2, 0.3])
#print "gradYlm(l=1, m=0, r=[0.0, 0.0, 0.3]): ", gradYlm(1, 0, [0.0, 0.0, 0.3])

def Xlm(l, m, pts):
  pts = numpy.asarray(pts)
  return numpy.cross(pts, gradYlm(l, m, pts)) / math.sqrt(l * (l + 1))


#print "Xlm(l=2, m=2, r=[0.1, 0.2, 0.3]): ", Xlm(2, 2, [0.1, 0.2, 0.3])
#print "Xlm(l=2, m=-2, r=[0.1, 0.2, 0.3]): ", Xlm(2, -2, [0.1, 0.2, 0.3])

def curlXlm(l, m, pts):
  pts = numpy.asarray(pts)
  fl = float(l)
  r, theta, phi = toSph(pts)

  res = numpy.zeros(r.shape + (3, ), dtype = 'complex')
  res += pts

  vals = fl * (fl + 1) * sph_harm(m, l, phi, theta) / (r * r)
  for i in range(3):
    res[...,i] *= vals

  res += gradYlm(l, m, pts)
  return -res / math.sqrt(l * (l + 1))

#print "curlXlm(l=2, m=2, r=[0.1, 0.2, 0.3]): ", curlXlm(2, 2, [0.1, 0.2, 0.3])

def aTM(l, m, k):
  kaIn = rtEpsIn * k * a
  kaOut = rtEpsOut * k * a
  kbOut = rtEpsOut * k * b
  res = jl(l, kaIn)
  return res / (jl(l, kaOut) - alt_jl(l, kbOut) * yl(l, kaOut) / alt_yl(l, kbOut))

#print "aTM_out(l=2, m=1, k=kTM021): ", aTM_out(2, 1, 4.3165735561686065791)

def bTM(l, m, k):
  kbOut = rtEpsOut * k * b
  return -alt_jl(l, kbOut) * aTM(l, m, k) / alt_yl(l, kbOut)

#print "bTM_out(l=2, m=1, k=kTM021): ", bTM_out(2, 1, 4.3165735561686065791)

def aTE(l, m, k):
  kaIn = rtEpsIn * k * a
  kaOut = rtEpsOut * k * a
  kbOut = rtEpsOut * k * b
  res = jl(l, kaIn)
  res /= jl(l, kaOut) - jl(l, kbOut) * yl(l, kaOut) / yl(l, kbOut)
  return res
  
def bTE(l, m, k):
  kbOut = rtEpsOut * k * b
  return -aTE(l, m, k) * jl(l, kbOut) / yl(l, kbOut)

#def bfieldTM(l, m, k, pt):
#  pt = numpy.asarray(pt, dtype = 'd')
#  r, theta, phi = list(toSph(pt))
#  if r < a:
#    return sph_jn(l, rtEpsIn * k * r)[0][-1] * (Xlm(l, m, pt) + (-1)**m * Xlm(l, -m, pt))
#  elif r < b:
#    res = aTM(l, m, k) * sph_jn(l, k * r)[0][-1]
#    res += bTM(l, m, k) * sph_yn(l, k * r)[0][-1]
#    res *= Xlm(l, m, pt) + (-1)**m * Xlm(l, -m, pt)
#    return res
#  else:
#    return numpy.zeros(3, dtype = 'd')

def efieldTM(l, m, k, pts, kind = "re"):
  pts = numpy.asarray(pts, dtype = 'd')
  r, theta, phi = toSph(pts)
  if kind == "re":
    realXlm = Xlm(l, m, pts) + (-1.0)**m * Xlm(l, -m, pts)
    realCurlXlm = curlXlm(l, m, pts) + (-1)**m * curlXlm(l, -m, pts)
  else:
    realXlm = -1.0j * (Xlm(l, m, pts) - (-1.0)**m * Xlm(l, -m, pts))
    realCurlXlm = -1.0j * (curlXlm(l, m, pts) - (-1)**m * curlXlm(l, -m, pts))

  res = numpy.zeros(r.shape + (3, ), dtype = 'd')

  indsIn = numpy.logical_and(r > 0.0, r < a)
  rIn = r[indsIn]
  rIn = rIn.reshape(rIn.shape + (1,))
  xIn = rtEpsIn * k * rIn
  res[indsIn] = rtEpsIn * djl(l, xIn) * numpy.cross(pts[indsIn], realXlm[indsIn]) / rIn
  res[indsIn] += jl(l, xIn) * realCurlXlm[indsIn] / k
  res[indsIn] /= epsIn

  indsOut = numpy.logical_and(r >= a, r < b)
  rOut = r[indsOut]
  rOut = rOut.reshape(rOut.shape + (1,))
  xOut = rtEpsOut * k * rOut
  rCrossXlm = numpy.cross(pts[indsOut], realXlm[indsOut])
  res[indsOut] = aTM(l, m, k) * (djl(l, xOut) * rCrossXlm / rOut + jl(l, xOut) * realCurlXlm[indsOut] / k)
  res[indsOut] += bTM(l, m, k) * (dyl(l, xOut) * rCrossXlm / rOut + yl(l, xOut) * realCurlXlm[indsOut] / k)

  return res


def efieldTE(l, m, k, pt, phase = "re"):
  pt = numpy.asarray(pt, dtype = 'd')
  r, theta, phi = toSph(pt)
  if phase == "re":
    realXlm = Xlm(l, m, pt) + (-1.0)**m * Xlm(l, -m, pt)
  else:
    realXlm = -1.0j * (Xlm(l, m, pt) - (-1.0)**m * Xlm(l, -m, pt))
  #print realXlm
  #print r, theta, phi

  res = numpy.zeros(r.shape + (3, ), dtype = 'd')

  indsIn = numpy.logical_and(r > 0.0, r < a)
  rIn = r[indsIn]
  rIn = rIn.reshape(rIn.shape + (1,))
  xIn = rtEpsIn * k * rIn
  res[indsIn] = jl(l, xIn) * realXlm[indsIn]

  indsOut = numpy.logical_and(r >= a, r < b)
  rOut = r[indsOut]
  rOut = rOut.reshape(rOut.shape + (1,))
  xOut = rtEpsOut * k * rOut
  res[indsOut] = (aTE(l, m, k) * jl(l, xOut) + bTE(l, m, k) * yl(l, xOut)) * realXlm[indsOut]

  return res

def efield(pol, n, l, m, phase, pt):
  k = findK(pol, n, l)
  print "k value: ", k
  if pol == "TE":
    #k = kTEs[n, l]
    return efieldTE(l, m, k, pt, phase)
  elif pol == "TM":
    #k = kTMs[n, l]
    return efieldTM(l, m, k, pt, phase)
  else:
    print "dsph_in_msph.efield(...): not a valid polarization."


def _octLowerBCsTMre(l, m):
  mOdd = m % 2
  lmmOdd = (l - m) % 2
  res = []
  # x = 0 bndry for real part
  if mOdd:
    res.append("pec")
  else:
    res.append("pmc")
  # y = 0 bndry for real part
  res.append("pmc")
  if lmmOdd:
    res.append("pec")
  else:
    res.append("pmc")
  return res

def _octLowerBCsTMim(l, m):
  res = _octLowerBCsTMre(l, m)
  # the x,y = 0 boundary conditions are opposite from the RealTM b.c. z = 0 is the same
  for i in range(2):
    if res[i] == "pmc":
      res[i] = "pec"
    else:
      res[i] = "pmc"
  return res

def _octLowerBCsTM(l, m, phase):
  if phase == "re":
    return _octLowerBCsTMre(l, m)
  elif phase == "im":
    return _octLowerBCsTMim(l, m)
  else:
    print "Not a valid phase, choose 're' or 'im'"

def _octLowerBCsTE(l, m ,phase):
  res = _octLowerBCsTM(l, m, phase)
  # all TE boundary conditions are the exact opposite from same-phase TM b.c.
  for i in range(3):
    if res[i] == "pmc":
      res[i] = "pec"
    else:
      res[i] = "pmc"
  return res

def octLowerBCs(pol, l, m, phase):
  """
  Return the -x,-y,-z boundary conditions given a spherical multipole solution
  to Maxwell's equations (assuming that the -x,-y,-z corner of the simulation is the
  origin.
  
    pol = 'TE' or 'TM'
    phase = 're' or 'im'
    l, m  = integers where l > 0 and |m| <= l (l = 0 are static solutions)
  """

  if pol == "TE":
    return _octLowerBCsTE(l, m, phase)
  elif pol == "TM":
    return _octLowerBCsTM(l, m, phase)
  else:
    print "Not a valid polarization. Should be 'TM' or 'TE'"


#testPt = [0.1, 0.2, 0.3]
#print numpy.cross(testPt, Xlm(2, 1, testPt) + (-1.0)**1 * Xlm(2, -1, testPt))

#print efieldTM(2, 1, 4.3165735561686065791, [0.1, 0.2, 0.3], "re")
#print bfieldTM(2, 1, 4.3165735561686065791, [0.1, 0.2, 0.3])
#print efieldTE(2, 2, 4.0593020005144315490, [0.1, 0.2, 0.3], "im")
#print efieldTE(2, 2, 4.0593020005144315490, [0.1, 0.2, 0.3], "re")



# low k-values from mathematica
# access via kTEs[n, l]
kTEs = numpy.zeros((3, 8), dtype = 'd')
# l = 1
kTEs[0, 1] = 3.0859803649032310262
kTEs[1, 1] = 5.5909053475831334718
kTEs[2, 1] = 8.1637117098749546912
# l = 2
kTEs[0, 2] = 4.0593020005144315490
kTEs[1, 2] = 6.7051972790870259596
kTEs[2, 2] = 9.3435078451641308394
# l = 3
kTEs[0, 3] = 5.0341266357880822580
kTEs[1, 3] = 7.8011739287676578998
kTEs[2, 3] = 10.500312224157252631
# l = 4
kTEs[0, 4] = 6.0114934306849569915
kTEs[1, 4] = 8.8834948012229946328
# l = 5
kTEs[0, 5] = 6.9894730426536395198
kTEs[1, 5] = 9.9542595209389566059
# l = 6
kTEs[0, 6] = 7.9658575396158286708
# l = 7
kTEs[0, 7] = 8.9389561607112238036

# access via kTMs[n, l]
kTMs = numpy.zeros((3, 8), dtype = 'd')
# l = 1
kTMs[0, 1] = 2.7914257502896397024
kTMs[1, 1] = 4.5810385029697019144
kTMs[2, 1] = 6.9275511065885668635
# l = 2
kTMs[0, 2] = 4.3165735561686065791
kTMs[1, 2] = 6.1047106593034911642
kTMs[2, 2] = 8.1812873773620704275
# l = 3
kTMs[0, 3] = 5.5881528173032377250
kTMs[1, 3] = 7.7462579349784372885
kTMs[2, 3] = 9.5156360470082189160
# l = 4
kTMs[0, 4] = 6.7149977850762035456
kTMs[1, 4] = 9.2802721914218748026
# l = 5
kTMs[0, 5] = 7.7719636368369787820
kTMs[1, 5] = 10.606364456388263491
# l = 6
kTMs[0, 6] = 8.7916393433331864203
# l = 7
kTMs[0, 7] = 9.7886040142934349323


fTMs = lightspeed * numpy.asarray(kTMs) / (2.0 * math.pi)
fTEs = lightspeed * numpy.asarray(kTEs) / (2.0 * math.pi)


def strToMode(str):
  pol = str.split("-")[0]
  n = int(str.split("-")[1])
  l = int(str.split("-")[2])
  m = int(str.split("-")[3])
  phase = str.split("-")[4]
  return pol, n, l, m, phase

def modeToVtk(resolution, modename, pathToH5Utils = ""):
  import tables

  delta = 0.5 / float(resolution)

  ex = numpy.zeros(3 * (resolution, ), dtype = 'd')
  ey = numpy.zeros(3 * (resolution, ), dtype = 'd')
  ez = numpy.zeros(3 * (resolution, ), dtype = 'd')

  pol, n, l, m, phase = strToMode(modename)
  print "Saving " + modename + " field"
  #print "Saving TM%d%d%d field" % (n, l, m)
  fx = "elecField_" + modename + "_x.h5"
  fy = "elecField_" + modename + "_y.h5"
  fz = "elecField_" + modename + "_z.h5"
  fvtk = "elecField_" + modename + ".vtk"
  h5x = tables.openFile(fx, 'w')
  h5y = tables.openFile(fy, 'w')
  h5z = tables.openFile(fz, 'w')

  pts = numpy.mgrid[0:0.5:delta, 0:0.5:delta, 0:0.5:delta]
  pts = numpy.rollaxis(pts, 0, pts.ndim)
  field = efield(pol, n, l, m, phase, pts)

  h5x.createArray("/", "data", field[...,0])
  h5y.createArray("/", "data", field[...,1])
  h5z.createArray("/", "data", field[...,2])
  h5x.close()
  h5y.close()
  h5z.close()
  # now use h5tovtk to get vtk file
  os.system(pathToH5Utils + "h5tovtk -o %s %s %s %s" % (fvtk, fx, fy, fz))
  #os.system("rm " + fx)
  #os.system("rm " + fy)
  #os.system("rm " + fz)


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



