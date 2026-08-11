package ffsqrt

import (
	"fmt"

	"github.com/ReneBoedker/algobra/finitefield/extfield"
)

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
