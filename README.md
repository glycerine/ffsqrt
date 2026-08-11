ffsqrt
======
based on https://github.com/ReneBoedker/algobra ,
to compute the square root of elements from a 
general finite field GF(p^k), where p is a prime number, 
and k an integer k>=1:

You don't need a different library. You don't need 
to reimplement anything field-specific. Tonelli–Shanks 
(the algorithm [libnum](https://github.com/hellman/libnum/blob/d90dc9ec5769bcadd483f98f0c71b587ceeb80f0/libnum/sqrtmod.py#L74) uses) 
only depends on F* being 
a cyclic group of order q-1, where q = p^k ; it never 
touches the internal representation of elements. 

So the same algorithm that finds square roots mod 
a prime also finds them in GF(p^k), as long as you 
have multiplication, inversion, and exponentiation, 
plus one element that generates F*. 

algobra's ff.Field/ff.Element interfaces give you those.

------
License: BSD 3-clause, the same as Go.
