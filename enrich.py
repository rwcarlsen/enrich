
import math
import matplotlib.pyplot as plt

kWhkg2GWdMT = 24000

def linspace(start = 0, end = 1, count = 1):
  vals = [0]*count
  step = float(end - start) / float(count)
  for i in range(count):
    vals[i] = start + i * step
  return vals

class Costs:
  def __init__(self, enrich = 133, convert = 6.5, mine = 100, dispose = 600):
    self.enrich = enrich
    self.convert = convert
    self.mine = mine
    self.dispose = dispose

class Matl:
  def __init__(self, mass = 1):
    self.Form = ''
    self.E = 0
    self.Enrich = 0
    self.Val = 0
    self.M = mass

  def ValPerKg(self):
    return self.Val / self.M

  def ValPerKWh(self):
    return self.Val / self.E

def mine(m, cost):
  m.Form = "U3O8"
  m.Enrich = 0.0071
  m.Val += cost * m.M
  return m

def convert(m, cost):
  m.Form = "U"
  UvTot = 3. * 238 / (3. * 238 + 8 * 16)
  m.M = UvTot * m.M
  m.Val += cost * m.M
  return m

def enrich(m, cost, tgt, tail):
  P, nswu = swu(m.M, m.Enrich, tgt, tail)

  m.Val += cost * nswu
  m.M = P
  m.Enrich = tgt

  return m

def burn(m, nBatches):
  # modify b1 formula to taste
  b1 = 1000 * m.Enrich

  # LRM discharge burnup
  bd = 2.0 * nBatches / (nBatches + 1) * b1

  m.E = bd * kWhkg2GWdMT * m.M
  return m

def dispose(m, cost):
  m.Val += cost * m.M
  return m

def swu(F, xf, xp, xt):
  P = (xf - xt) / (xp - xt) * F
  T = F - P
  return P, P * potential(xp) + T * potential(xt) - F * potential(xf)

def potential(x):
  return (1 - 2 * x) * math.log((1 - x) / x)

def onceThrough(m, costs = None, prod = 0.035, tail = 0.003, nBatches = 3):
  if costs is None:
    costs = Costs()

  mine(m, costs.mine)
  convert(m, costs.convert)
  enrich(m, costs.enrich, prod, tail)
  burn(m, nBatches)
  dispose(m, costs.dispose)

def frontend(m, costs = None, prod = 0.035, tail = 0.003, nBatches = 3):
  if costs is None:
    costs = Costs()

  mine(m, costs.mine)
  convert(m, costs.convert)
  enrich(m, costs.enrich, prod, tail)
  burn(m, nBatches)

def backend(m, costs = None, prod = 0.035, tail = 0.003, nBatches = 3):
  if costs is None:
    costs = Costs()

  burn(m, nBatches)
  dispose(m, costs.dispose)

def vary_dispose():
  maxenr = .2
  natenr = .0071
  enrichments = linspace(natenr, maxenr, 10000)
  disp_costs = linspace(0, 1000, 10)
  for cost in disp_costs:
    valkg = []
    valkwh = []
    for e in enrichments:
      m = Matl()
      onceThrough(m, Costs(dispose = cost), prod = e)

      valkg.append(m.ValPerKg())
      valkwh.append(m.ValPerKWh())

    plt.plot(enrichments, valkwh, label = str(cost))
  plt.show()

def vary_mining():
  maxenr = .2
  natenr = .0071
  enrichments = linspace(natenr, maxenr, 10000)
  mine_costs = linspace(0, 200, 10)
  for cost in mine_costs:
    valkg = []
    valkwh = []
    for e in enrichments:
      m = Matl()
      onceThrough(m, Costs(mine = cost), prod = e)

      valkg.append(m.ValPerKg())
      valkwh.append(m.ValPerKWh())

    plt.plot(enrichments, valkwh, label = str(cost))
  plt.show()


def swuplot():
  maxenr = .2
  natenr = .0071
  enrichments = linspace(natenr, maxenr, 10000)

  swus = []

  for e in enrichments:
    p, s = swu(1, natenr, e, 0.003)
    # s / p for swus/kg-product or just s for swus/kg-feed
    swus.append(s / p)

  plt.plot(enrichments, swus)
  plt.show()

def only_dispose():
  maxenr = .1
  natenr = .0071
  enrichments = linspace(natenr, maxenr, 10000)

  kgcosts = []
  kwhcosts = []

  for e in enrichments:
    m = Matl()
    backend(m, prod = e)
    # s / p for swus/kg-product or just s for swus/kg-feed
    kgcosts.append(m.ValPerKg())
    kwhcosts.append(m.ValPerKWh())

  plt.plot(enrichments, kwhcosts)
  plt.show()

def only_enrich():
  maxenr = .1
  natenr = .0071
  enrichments = linspace(natenr, maxenr, 10000)

  kgcosts = []
  kwhcosts = []

  for e in enrichments:
    m = Matl()
    frontend(m, prod = e)
    # s / p for swus/kg-product or just s for swus/kg-feed
    kgcosts.append(m.ValPerKg())
    kwhcosts.append(m.ValPerKWh())

  plt.plot(enrichments, kwhcosts)
  plt.show()

if __name__ == '__main__':
  vary_dispose()
  vary_mining()
  #only_enrich()
  #swuplot()

