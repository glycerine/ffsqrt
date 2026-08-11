ffsqrt
======

Was based on https://github.com/ReneBoedker/algobra ,
but that has a very bad bug in it (see https://github.com/ReneBoedker/algobra/issues/2 ), so, for now, based on my fork https://github.com/glycerine/algobra :

to compute the square root of elements from a 
general finite field GF(p^k), where p is a prime number, 
and k an integer k>=1:

You don't need a different library; algobra (fixed) suffices.
Tonelli–Shanks (the algorithm [libnum](https://github.com/hellman/libnum/blob/d90dc9ec5769bcadd483f98f0c71b587ceeb80f0/libnum/sqrtmod.py#L74) uses) 
only depends on F* being a cyclic group of 
order q-1, where q = p^k ; it never 
touches the internal representation of elements. 

So the same algorithm that finds square roots mod 
a prime also finds them in GF(p^k), as long as you 
have multiplication, inversion, and exponentiation, 
plus a quadratic non-residue to seed the algorithm. 

algobra's ff.Field/ff.Element interfaces give you those,
so we just use them here.

Two things the enclosed sqrt.go does better, over a naive port of libnum:

* Char-2 fields are trivial. In GF(2^k) every element is a square (squaring is the Frobenius automorphism), and the inverse is just a^(q/2). No Tonelli–Shanks loop is needed.

* No random search for a non-residue. libnum finds a non-residue by testing random elements. Here, we first try f.MultGenerator() and verify g^((q-1)/2) != 1. If that helper returns a non-generator, we fall back to a deterministic search for a usable non-residue.


------
License: BSD 3-clause, the same as Go.
