
import math
import matplotlib.pyplot as plt

kWhkg2GWdMT = 24000

class Matl:
  def __init__(self):
    self.Form = ''
    self.E = 0
    self.Enrich = 0
    self.Val = 0
    self.M = 0

  def ValPerKg(self):
    return self.Val / self.M

  def ValPerKWh(self):
    return self.Val / self.E

def mine(m):
  m.Form = "U3O8"
  m.Enrich = 0.0071
  m.Val += 100 * m.M
  return m

def convert(m):
  m.Form = "U"
  UvTot = 3. * 238 / (3. * 238 + 8 * 16)
  m.M = UvTot * m.M
  m.Val += 6.5 * m.M
  return m

def enrich(m, tgt, tail):
  P, nswu = swu(m.M, m.Enrich, tgt, tail)

  m.Val += 133 * nswu
  m.Form = "U"
  m.M = P
  m.Enrich = tgt

  return m

def burn(m, nBatches):
  # modify b1 formula to taste
  b1 = 1000 * m.Enrich

  # LRM discharge burnup
  n = nBatches
  bd = 2 * n / (n + 1) * b1

  m.E = bd * kWhkg2GWdMT * m.M
  return m

def dispose(m):
  m.Val += 2000 * m.M
  return m

def swu(F, xf, xp, xt):
  P = (xf - xt) / (xp - xt) * F
  T = F - P
  return P, P * potential(xp) + T * potential(xt) - F * potential(xf)

def potential(x):
  return (1 - 2 * x) * math.log((1 - x) / x)

def onceThrough(m, e):
  mine(m)
  convert(m)
  enrich(m, e, 0.003)
  burn(m, 3)
  dispose(m)

if __name__ == '__main__':

  maxenr = .2
  natenr = .0071
  delta = (maxenr - natenr) / 100

  enrichments = []
  valkg = []
  valkwh = []

  e = natenr
  while e <= maxenr:
    e += delta

    m = Matl()
    m.M = 1
    onceThrough(m, e)

    enrichments.append(m.Enrich)
    valkg.append(m.ValPerKg())
    valkwh.append(m.ValPerKWh())

  plt.plot(enrichments, valkg)

