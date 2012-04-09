
package main

import (
  "fmt"
  "math"
)

const (
  kWhkg2GWdMT = 24000
)

func main() {

  maxenr := .2
  for e := 0.0071; e <= maxenr; e += (maxenr - .0071) / 10 {
    m := &Matl{M:1}
    onceThrough(m, Enrichment(e))
    fmt.Println("\ncost:")
    fmt.Println(m)
    fmt.Println(m.ValPerKWh())
  }

}

func onceThrough(m *Matl, e Enrichment) {
  mine(m)
  convert(m)
  enrich(m, e, 0.003)
  burn(m, 3)
  dispose(m)
}

type (
  Chemform string
  Enrichment float64 // mass frac 235
  Burnup float64 // GWD/MTU
  Value float64 // $
  Energy float64 // kWh
  Mass float64 // kg
)

type Matl struct {
  Form Chemform
  E Energy
  Enrich Enrichment
  Val Value
  M Mass
}

func (m *Matl) ValPerKg() Value {
  return m.Val / Value(m.M)
}

func (m *Matl) ValPerKWh() Value {
  return m.Val / Value(m.E)
}

func mine(m *Matl) *Matl {
  m.Form = "U3O8"
  m.Enrich = 0.0071
  m.Val += 100 * Value(m.M)
  return m
}

func convert(m *Matl) *Matl {
  m.Form = "U"
  UvTot := 3. * 238 / (3. * 238 + 8 * 16)
  m.M = Mass(UvTot) * m.M
  m.Val += 6.5 * Value(m.M)
  return m
}

func enrich(m *Matl, tgt, tail Enrichment) *Matl {
  P, nswu := swu(m.M, m.Enrich, tgt, tail)

  m.Val += 133 * Value(nswu)
  m.Form = "U"
  m.M = P
  m.Enrich = tgt

  return m
}

func burn(m *Matl, nBatches int) *Matl {
  // modify b1 formula to taste
  b1 := Burnup(1000) * Burnup(m.Enrich)

  // LRM discharge burnup
  n := Burnup(nBatches)
  bd := 2 * n / (n + 1) * b1

  m.E = Energy(bd) * kWhkg2GWdMT * Energy(m.M)
  return m
}

func dispose(m *Matl) *Matl {
  m.Val += 2000 * Value(m.M)
  return m
}

func swu(F Mass, xf, xp, xt Enrichment) (Mass, float64) {
  P := Mass((xf - xt) / (xp - xt)) * F
  T := F - P
  return P, float64(P) * potential(xp) + float64(T) * potential(xt) - float64(F) * potential(xf)
}

func potential(enrich Enrichment) float64 {
  x := float64(enrich)
  return (1 - 2 * x) * math.Log((1 - x) / x)
}

