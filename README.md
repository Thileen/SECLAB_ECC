# SECLAB_ECC

SECLAB_ECC is a collection with some 32-bit scalar multiplication over a finite field including countermeasures for the elliptic curve in Weierstrass and Twisted Edwards form. The scalar multiplication is implemented with the binary method. 

These following countermeasures are implemented:

* Double-And-Add-Always 
* Montgomery Ladder
* Scalar Randomization
* Randomized Scalar Splitting
* Random Projective Coordinates
* Randomized Sliding Window

The following domain parameter can be used for the elliptic curve in Weierstrass form:

* p_192
* p_224
* p_256
* p_384
* p_521
* secp128r1
* secp192k1
* secp256k1
* secp256r1
* brainpoolP256r1
* brainpoolP320r1
* brainpoolP384r1
* brainpoolP512r1

and for the elliptic cuve in Twisted Edwards form 

* ed25519

The reduction is calculated by the Barrett-Reduction, only for the domain parameter p_256 exist the fast reduction mathod.

The testbench.c provide a cup of examples how to use the implementations

# ToDo:
* implement more domain parameter
* implement ECC scheme (ECDH, ECDIAS, ECDSA)
