// Package ffsqrt implements modular square roots over an arbitrary finite
// field GF(p^k), using a generalization of the Tonelli-Shanks algorithm.
//
// Tonelli-Shanks only depends on the multiplicative group F* being cyclic of
// order q-1 (q = p^k) — it never looks at the internal representation of
// field elements. That means the exact same algorithm that computes square
// roots mod a prime p also computes square roots in GF(p^k): you just need
// field multiplication, inversion, and exponentiation, plus a quadratic
// non-residue to seed the algorithm. A multiplicative generator is guaranteed
// to be one, but this package verifies that seed and falls back to a
// deterministic search if needed.
//
// Built against github.com/glycerine/algobra v0.2.2-jea's ff.Field / ff.Element
// interfaces (github.com/glycerine/algobra/finitefield/ff), which are
// implemented by extfield.Field (arbitrary p^k), primefield.Field (GF(p)),
// and binfield.Field (GF(2^k)) alike, so this code works unmodified against
// any of them.
package ffsqrt

import (
	"errors"

	"github.com/glycerine/algobra/finitefield/ff"
)

// ErrNonResidue is returned when a has no square root in the field, i.e. a is
// a quadratic non-residue.
var ErrNonResidue = errors.New("ffsqrt: element is not a quadratic residue")

var errNoNonResidue = errors.New("ffsqrt: field has no quadratic non-residue")

// Sqrt returns a square root of a in its field f. If a is not a square, it
// returns ErrNonResidue. The other root (if any) is -root.
//
// f must be the field that a is defined over.
func Sqrt(f ff.Field, a ff.Element) (ff.Element, error) {
	if a.IsZero() {
		return f.Zero(), nil
	}

	q := f.Card() // q = p^k

	// Characteristic 2: squaring is the Frobenius endomorphism x -> x^p,
	// which is a field automorphism, so *every* element is a square, and
	// it's cheap to invert directly. Since x^q = x for all x (Fermat),
	// squaring applied k times is the identity, so squaring applied k-1
	// more times inverts a single squaring: sqrt(a) = a^(q/2).
	if f.Char() == 2 {
		return a.Pow(q / 2), nil
	}

	qm1 := q - 1
	half := qm1 / 2 // (q-1)/2

	// Euler's criterion: a is a nonzero square iff a^((q-1)/2) == 1.
	if euler := a.Pow(half); !euler.IsOne() {
		return nil, ErrNonResidue
	}

	// Factor q-1 = Q * 2^S with Q odd.
	Q := qm1
	S := uint(0)
	for Q%2 == 0 {
		Q /= 2
		S++
	}

	// Fast path: q ≡ 3 (mod 4), i.e. S == 1. Then sqrt(a) = a^((q+1)/4),
	// which in terms of Q (since q-1 = 2Q here) is a^((Q+1)/2).
	if S == 1 {
		return a.Pow((Q + 1) / 2), nil
	}

	// General Tonelli-Shanks. We need one quadratic non-residue to seed
	// the algorithm. A generator of F* is guaranteed to be one, but some
	// ff.Field implementations have returned non-generators here, so verify
	// the seed before trusting it.
	z, err := findNonResidue(f, half)
	if err != nil {
		return nil, err
	}

	c := z.Pow(Q)
	t := a.Pow(Q)
	R := a.Pow((Q + 1) / 2)
	M := S

	for {
		if t.IsOne() {
			return R, nil
		}

		// Find the least i, 0 < i < M, with t^(2^i) == 1.
		i := uint(0)
		temp := t
		for !temp.IsOne() {
			temp = temp.Times(temp) // Times is non-mutating (unlike Mult)
			i++
			if i == M {
				// Shouldn't happen: Euler's criterion already
				// guaranteed a is a residue.
				return nil, ErrNonResidue
			}
		}

		b := c
		for j := uint(0); j < M-i-1; j++ {
			b = b.Times(b)
		}

		M = i
		c = b.Times(b)
		t = t.Times(c)
		R = R.Times(b)
	}
}

func findNonResidue(f ff.Field, half uint) (ff.Element, error) {
	if z := f.MultGenerator(); isNonResidue(z, half) {
		return z, nil
	}

	// In prime fields this usually finds a non-residue without allocating all
	// field elements. In extension fields of even degree the prime subfield may
	// contain only squares, so this is only a cheap first fallback.
	for n := uint(2); n < f.Char(); n++ {
		if z := f.ElementFromUnsigned(n); isNonResidue(z, half) {
			return z, nil
		}
	}

	for _, z := range f.Elements() {
		if isNonResidue(z, half) {
			return z, nil
		}
	}

	return nil, errNoNonResidue
}

func isNonResidue(a ff.Element, half uint) bool {
	if a == nil || a.IsZero() || a.Err() != nil {
		return false
	}

	euler := a.Pow(half)
	return euler.Err() == nil && !euler.IsOne()
}
