
package main

import (
  "fmt"
  "math"
)

func main() {
  m := &Matl{M:1}

  mine(m)
  convert(m)
  enrich(m, 0.035, 0.003)

  fmt.Println(m)
  fmt.Println(m.ValPerKg())
}

type (
  Chemform string
  Enrichment float64
  Value float64
  Mass float64
)

type Matl struct {
  Form Chemform
  Enrich Enrichment // mass frac 235
  Val Value // $/kg
  M Mass // kg
}

func (m *Matl) ValPerKg() Value {
  return m.Val / Value(m.M)
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


