package ffsqrt

import (
	"fmt"
	"testing"

	"github.com/glycerine/algobra/finitefield/primefield"
)

var _ = fmt.Printf

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

// TestSqrt_ExhaustiveSweep_AllPrimesUnderBound exhaustively checks every
// element of GF(p) for every prime p below sweepBound: for each a, Sqrt(a^2)
// must return some root whose square is a^2. This exists because the
// original hand-picked prime list (7, 11, 13, 17, 19, 29, 31, 97, ...) missed
// a real bug (algobra v0.2.1's primefield.Field.MultGenerator occasionally
// returning a non-generator, e.g. for GF(17), GF(97), GF(65537)) — this test
// sweeps broadly instead of trusting a short list of hand-picked primes.
//
// It's the most expensive test in the suite (a few thousand field
// constructions, each exhaustively iterated), so it's skipped under
// `go test -short`.
func TestSqrt_ExhaustiveSweep_AllPrimesUnderBound(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping exhaustive prime sweep in -short mode")
	}

	const sweepBound = 3000
	for p := uint(2); p < sweepBound; p++ {
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
