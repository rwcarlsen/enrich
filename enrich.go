
package main

import (
  "fmt"
  "math"
)

const (
  kWhkg2GWdMT = 24000
)

func main() {
  m := &Matl{M:1}

  mine(m)
  convert(m)
  enrich(m, 0.035, 0.003)
  dispose(m)

  fmt.Println(m)
  fmt.Println(m.ValPerKg())
  fmt.Println(m.ValPerKWh(bu))
}

func onceThrough(m *Matl) {
  mine(m)
  convert(m)
  enrich(m, 0.035, 0.003)
  burn(m, 3)
  dispose(m)

  fmt.Println(m)
  fmt.Println(m.ValPerKg())
  fmt.Println(m.ValPerKWh())
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

func (m *Matl) ValPerKWh(bu Burnup) Value {
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
  xf := m.Enrich
  xp := tgt
  xt := tail
  F := m.M
  P := Mass((xf - xt) / (xp - xt)) * F
  T := F - P

  swu := float64(P) * potential(xp) + float64(T) * potential(xt) - float64(F) * potential(xf)

  m.Val += 133 * Value(swu)
  m.Form = "U"
  m.M = P
  m.Enrich = xp

  return m
}

func potential(enrich Enrichment) float64 {
  x := float64(enrich)
  return (1 - 2 * x) * math.Log((1 - x) / x)
}

func dispose(m *Matl) *Matl {
  m.Val += 2000 * Value(m.M)
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

