package ffsqrt

import (
	"testing"

	"github.com/glycerine/algobra/finitefield/primefield"
)

const primeSweepBound = 3000

// isPrime is a simple trial-division primality test. It only needs to be
// correct, not fast, since it's only used to build the list of primes to
// sweep over below.
func isPrime(n uint) bool {
	if n < 2 {
		return false
	}
	for d := uint(2); d*d <= n; d++ {
		if n%d == 0 {
			return false
		}
	}
	return true
}

func distinctPrimeFactors(n uint) []uint {
	var factors []uint
	for d := uint(2); d*d <= n; d++ {
		if n%d != 0 {
			continue
		}
		factors = append(factors, d)
		for n%d == 0 {
			n /= d
		}
	}
	if n > 1 {
		factors = append(factors, n)
	}
	return factors
}

// TestPrimeFieldMultGenerator_Sweep_AllPrimesUnderBound checks the exact
// contract Sqrt used to rely on: primefield.Field.MultGenerator must return a
// generator of GF(p)*. A nonzero element g is a generator iff g^(p-1)=1 and
// g^((p-1)/r)!=1 for every prime factor r of p-1.
func TestPrimeFieldMultGenerator_Sweep_AllPrimesUnderBound(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping exhaustive prime generator sweep in -short mode")
	}

	for p := uint(2); p < primeSweepBound; p++ {
		if !isPrime(p) {
			continue
		}

		field, err := primefield.Define(p)
		if err != nil {
			t.Fatalf("primefield.Define(%d): %v", p, err)
		}

		g := field.MultGenerator()
		if g.IsZero() {
			t.Errorf("p=%d: MultGenerator() = 0, want generator of GF(p)*", p)
			continue
		}
		if got := g.Pow(p - 1); !got.IsOne() {
			t.Errorf("p=%d: generator candidate %s has g^(p-1) = %s, want 1", p, g, got)
		}

		for _, factor := range distinctPrimeFactors(p - 1) {
			if got := g.Pow((p - 1) / factor); got.IsOne() {
				t.Errorf("p=%d: generator candidate %s has g^((p-1)/%d) = 1", p, g, factor)
			}
		}
	}
}

// TestSqrt_ExhaustiveSweep_AllPrimesUnderBound exhaustively checks every
// element of GF(p) for every prime p below sweepBound: for each a, Sqrt(a^2)
// must return some root whose square is a^2. This complements the direct
// MultGenerator sweep above by checking Sqrt's full public behavior over many
// prime fields instead of trusting a short list of hand-picked primes.
//
// It's the most expensive test in the suite (a few thousand field
// constructions, each exhaustively iterated), so it's skipped under
// `go test -short`.
func TestSqrt_ExhaustiveSweep_AllPrimesUnderBound(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping exhaustive prime sweep in -short mode")
	}

	for p := uint(2); p < primeSweepBound; p++ {
		if !isPrime(p) {
			continue
		}
		field, err := primefield.Define(p)
		if err != nil {
			t.Fatalf("primefield.Define(%d): %v", p, err)
		}

		for _, a := range field.Elements() {
			square := a.Times(a)

			root, err := Sqrt(field, square)
			if err != nil {
				t.Errorf("p=%d: Sqrt(%s) [known square of %s] returned unexpected error: %v", p, square, a, err)
				continue
			}
			if got := root.Times(root); !got.Equal(square) {
				t.Errorf("p=%d: Sqrt(%s) = %s, but %s^2 = %s (want %s)", p, square, root, root, got, square)
			}
		}
	}
}
