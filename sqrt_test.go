package ffsqrt

import (
	"errors"
	"fmt"
	"testing"

	"github.com/ReneBoedker/algobra/finitefield/binfield"
	"github.com/ReneBoedker/algobra/finitefield/extfield"
	"github.com/ReneBoedker/algobra/finitefield/ff"
	"github.com/ReneBoedker/algobra/finitefield/primefield"
)

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

// assertSquaresRoundTrip checks, for every element a in the field, that
// Sqrt(a^2) returns some root r with r^2 == a^2. It does not require
// r == a, since square roots are only unique up to sign.
func assertSquaresRoundTrip(t *testing.T, field ff.Field) {
	t.Helper()

	for _, a := range field.Elements() {
		want := a.Times(a)

		root, err := Sqrt(field, want)
		if err != nil {
			t.Errorf("Sqrt(%s) [known square of %s] returned unexpected error: %v", want, a, err)
			continue
		}

		if got := root.Times(root); !got.Equal(want) {
			t.Errorf("Sqrt(%s) = %s, but %s^2 = %s (want %s)", want, root, root, got, want)
		}
	}
}

// assertNonResiduesRejected checks that Sqrt succeeds exactly on the values
// that are actually squares in the field, and fails with ErrNonResidue on
// everything else.
func assertNonResiduesRejected(t *testing.T, field ff.Field) {
	t.Helper()

	isSquare := make(map[string]bool)
	for _, a := range field.Elements() {
		isSquare[a.Times(a).String()] = true
	}

	for _, a := range field.Elements() {
		_, err := Sqrt(field, a)
		square := isSquare[a.String()]

		switch {
		case square && err != nil:
			t.Errorf("Sqrt(%s): got error %v, but %s IS a square", a, err, a)
		case !square && err == nil:
			t.Errorf("Sqrt(%s): got no error, but %s is NOT a square", a, a)
		case !square && !errors.Is(err, ErrNonResidue):
			t.Errorf("Sqrt(%s): got error %v, want ErrNonResidue", a, err)
		}
	}
}

// assertBothRootsDistinct checks that for a nonzero square in an odd
// characteristic field, r and -r are two distinct square roots (the only
// two, though this helper doesn't prove uniqueness by itself).
func assertBothRootsDistinct(t *testing.T, field ff.Field) {
	t.Helper()

	if field.Char() == 2 {
		return // -r == r in characteristic 2, nothing to check
	}

	for _, a := range field.Elements() {
		if a.IsZero() {
			continue
		}
		square := a.Times(a)

		root, err := Sqrt(field, square)
		if err != nil {
			t.Fatalf("Sqrt(%s) returned unexpected error: %v", square, err)
		}

		neg := root.Neg()
		if root.Equal(neg) {
			t.Errorf("root %s of nonzero square %s equals its own negation", root, square)
		}
		if got := neg.Times(neg); !got.Equal(square) {
			t.Errorf("-root (%s) does not square back to %s (got %s)", neg, square, got)
		}
	}
}

// runExhaustiveSuite runs the full set of exhaustive checks against field.
func runExhaustiveSuite(t *testing.T, name string, field ff.Field) {
	t.Run(name, func(t *testing.T) {
		t.Run("SquaresRoundTrip", func(t *testing.T) {
			assertSquaresRoundTrip(t, field)
		})
		t.Run("NonResiduesRejected", func(t *testing.T) {
			assertNonResiduesRejected(t, field)
		})
		t.Run("BothRootsDistinct", func(t *testing.T) {
			assertBothRootsDistinct(t, field)
		})
	})
}

// ---------------------------------------------------------------------
// Zero element
// ---------------------------------------------------------------------

func TestSqrt_Zero(t *testing.T) {
	for _, card := range []uint{2, 3, 5, 9, 16} {
		field, err := extfield.Define(card)
		if err != nil {
			t.Fatalf("extfield.Define(%d): %v", card, err)
		}

		root, err := Sqrt(field, field.Zero())
		if err != nil {
			t.Fatalf("Sqrt(0) in GF(%d) returned error: %v", card, err)
		}
		if !root.IsZero() {
			t.Errorf("Sqrt(0) in GF(%d) = %s, want 0", card, root)
		}
	}
}

// ---------------------------------------------------------------------
// GF(p): exhaustive, both p ≡ 1 (mod 4) and p ≡ 3 (mod 4), plus p = 2.
//
// p ≡ 3 (mod 4) exercises the S == 1 fast path.
// p ≡ 1 (mod 4) exercises the general Tonelli-Shanks loop (S > 1); 97-1 =
// 32*3, giving S = 5, which forces several iterations of the outer loop.
// ---------------------------------------------------------------------

func TestSqrt_PrimeFields_Exhaustive(t *testing.T) {
	primes := []uint{
		2,  // characteristic 2 special case
		3,  // 3 mod 4 == 3
		5,  // 5 mod 4 == 1
		7,  // 7 mod 4 == 3
		11, // 11 mod 4 == 3
		13, // 13 mod 4 == 1
		17, // 17 mod 4 == 1 (16 = 1*2^4, S=4)
		19, // 19 mod 4 == 3
		29, // 29 mod 4 == 1
		31, // 31 mod 4 == 3
		97, // 97 mod 4 == 1 (96 = 3*2^5, S=5)
	}

	for _, p := range primes {
		field, err := primefield.Define(p)
		if err != nil {
			t.Fatalf("primefield.Define(%d): %v", p, err)
		}
		runExhaustiveSuite(t, "GF_"+field.String(), field)
	}
}

// ---------------------------------------------------------------------
// GF(p^k), k > 1: exhaustive over cardinalities also used by algobra's own
// test suite, so we know Conway polynomials are available for them.
// ---------------------------------------------------------------------

func TestSqrt_ExtensionFields_Exhaustive(t *testing.T) {
	cardinalities := []uint{4, 8, 9, 16, 25, 49, 64}

	for _, q := range cardinalities {
		field, err := extfield.Define(q)
		if err != nil {
			t.Fatalf("extfield.Define(%d): %v", q, err)
		}
		runExhaustiveSuite(t, "GF_"+field.String(), field)
	}
}

// ---------------------------------------------------------------------
// GF(2^k) via the dedicated binfield package (char 2 fast path, exercised
// through a different underlying implementation than extfield).
// ---------------------------------------------------------------------

func TestSqrt_BinaryFields_Exhaustive(t *testing.T) {
	cardinalities := []uint{4, 8, 16, 32, 64, 256}

	for _, q := range cardinalities {
		field, err := binfield.Define(q)
		if err != nil {
			t.Fatalf("binfield.Define(%d): %v", q, err)
		}
		runExhaustiveSuite(t, "GF_"+field.String(), field)
	}
}

// ---------------------------------------------------------------------
// Known hand-checked values, as a readable sanity check independent of the
// exhaustive property tests above.
//
// In GF(7): 1^2=1, 2^2=4, 3^2=2, 4^2=2, 5^2=4, 6^2=1.
// QRs = {1, 2, 4}; non-residues = {3, 5, 6}.
// ---------------------------------------------------------------------

func TestSqrt_GF7_KnownValues(t *testing.T) {
	field, err := primefield.Define(7)
	if err != nil {
		t.Fatalf("primefield.Define(7): %v", err)
	}

	elem := func(v int) ff.Element {
		e, err := field.Element(v)
		if err != nil {
			t.Fatalf("field.Element(%d): %v", v, err)
		}
		return e
	}

	residues := map[int][2]int{
		1: {1, 6},
		2: {3, 4},
		4: {2, 5},
	}
	for a, roots := range residues {
		root, err := Sqrt(field, elem(a))
		if err != nil {
			t.Errorf("Sqrt(%d) in GF(7) returned error: %v", a, err)
			continue
		}
		got := 0
		for v := 0; v < 7; v++ {
			if root.Equal(elem(v)) {
				got = v
				break
			}
		}
		if got != roots[0] && got != roots[1] {
			t.Errorf("Sqrt(%d) in GF(7) = %d, want %d or %d", a, got, roots[0], roots[1])
		}
	}

	for _, a := range []int{3, 5, 6} {
		_, err := Sqrt(field, elem(a))
		if !errors.Is(err, ErrNonResidue) {
			t.Errorf("Sqrt(%d) in GF(7): got err=%v, want ErrNonResidue", a, err)
		}
	}
}

// ---------------------------------------------------------------------
// Randomized property test over a larger field, where full exhaustion
// would be needlessly slow. Complements the exhaustive suites above.
// ---------------------------------------------------------------------

func TestSqrt_RandomizedProperty_LargerField(t *testing.T) {
	field, err := primefield.Define(65537) // large prime, 65536 = 2^16 * 1 => S=16
	if err != nil {
		t.Fatalf("primefield.Define(65537): %v", err)
	}

	const trials = 500
	for i := 0; i < trials; i++ {
		a := field.RandElement()
		square := a.Times(a)

		root, err := Sqrt(field, square)
		if err != nil {
			t.Fatalf("trial %d: Sqrt(%s) [square of %s] returned error: %v", i, square, a, err)
		}
		if got := root.Times(root); !got.Equal(square) {
			t.Fatalf("trial %d: Sqrt(%s) = %s, but %s^2 = %s", i, square, root, root, got)
		}
	}
}

// ---------------------------------------------------------------------
// Field mismatch / defensive behavior: Sqrt should not panic when handed
// elements, and should behave consistently across repeated calls with the
// same input (determinism), even though internally MultGenerator/RandElement
// may be involved on other code paths.
// ---------------------------------------------------------------------

func TestSqrt_Deterministic(t *testing.T) {
	field, err := primefield.Define(97)
	if err != nil {
		t.Fatalf("primefield.Define(97): %v", err)
	}

	a, err := field.Element(10)
	if err != nil {
		t.Fatalf("field.Element(10): %v", err)
	}
	square := a.Times(a)

	first, err := Sqrt(field, square)
	if err != nil {
		t.Fatalf("Sqrt returned error: %v", err)
	}
	for i := 0; i < 20; i++ {
		again, err := Sqrt(field, square)
		if err != nil {
			t.Fatalf("call %d: Sqrt returned error: %v", i, err)
		}
		if !again.Equal(first) {
			t.Errorf("call %d: Sqrt(%s) = %s, want %s (same as first call)", i, square, again, first)
		}
	}
}

func Example() {
	// extfield.Define handles any prime power p^k, including plain
	// primes (k=1). Here: GF(3^4) = GF(81).
	field, err := extfield.Define(81)
	if err != nil {
		panic(err)
	}

	a, err := field.Element([]int{1, 2, 0, 1}) // a = 1 + 2x + 0x^2 + x^3
	if err != nil {
		panic(err)
	}

	square := a.Times(a) // guaranteed to be a quadratic residue

	root, err := Sqrt(field, square)
	if err != nil {
		panic(err) // shouldn't happen: square is a square by construction
	}

	fmt.Println(root.Times(root).Equal(square))
	// Output: true
}
