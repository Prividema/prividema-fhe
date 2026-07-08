# OnionPIR

[OnionPIR](https://eprint.iacr.org/2021/1081) is a private information retrieval (PIR) protocol.
As such, it allows a server to process encrypted queries for access
to specific database positions, without the server being able to
gain any insight as to the positions being accessed.

The prividema-fhe library provides an example implementation of OnionPIR,
based for the most part on the implementation presented as an example in
the [Hypets 2024 tutorials](https://github.com/tfhe/hypets2024_tutorial).
The implementation in said repository is a loose adaptation of the [original paper](https://eprint.iacr.org/2021/1081.pdf)
with modifications that leverage the performance improvements that the bivariate
decomposition permits.

## Packed GLWEGadgets {#packed_glwegadget}

The original OnionPIR protocol uses dot products of BFV ciphertexts with plaintext database entries.
The adapted version used in the Hypets tutorial uses bivariate representation for the database data,
and in that line uses the half-external product to achieve the same goal.

In order to lower the amount of data that needs to be transmitted to the server,
the OnionPIR protocol uses what it calls "query packing".
Using it, all the different BFV ciphertexts that define the query are _packed_ into
a single ciphertext.
The fact that allows this packing to be performed is that in OnionPIR,
the query ciphertexts have binary values (0 or 1) yet are represented as polynomials of high
degree (usually 1024 or more).
Thus, in the BFV ciphertexts that the protocol requires (for the products),
all but one of the coefficients in the polynomials remain zeroed and unused.
Then, for data transmission between the client and server, one can use a single polynomial's
coefficients to _pack_ all the coefficients that describe all the BFV ciphertexts that
need to be used, and encrypt it in a GLWE as usual.
The server can later on use an operation called "homomorphic trace" or "trace expansion" (depending on the reference)
to homomorphically expand a single GLWE ciphertext into a collection of GLWE ciphertexts that
each encrypt the

More formally, one can say that trace expansion is a process where a GLWE ciphertext of \f$ a(x) = a_0 + a_1 \cdot x + ... + a_d \cdot x^d\f$ is converted into a set of GLWE ciphertexts of \f$ a_0, a_1, \dots, a_d\f$ (possibly multiplied by a constant that depends on \f$d\f$).

Given that the half-external product uses a GLWEGadget, however, this is not directly usable.
For that reason, a way to pack GLWEGadgets into GLWE ciphertexts becomes necessary.
Leveraging the fact that GLWEGadgets are nothing more than a collection of GLWEs themeselves,
one can reuse the trace expansion algorithm if the original polynomial to encrypt is as follows:

Let the plaintext be \f$ a(x) = a_0 + ... + a_d \cdot x^d \f$, \f$N\f$ be the polynomial degree and \f$\tilde{\ell}\f$ the GLWEGadget depth parameter.
The GLWEGadget can be packed in a GLWE if and only if \f$ d \cdot \tilde{\ell} \le N\f$, in which case, its packed version is a(ny) encryption of \f$a_p(x) = A_1(x) + A_2(x) \cdot x^d + ... + A_{\tilde{\ell}}  \cdot  x^{(\tilde{\ell}-1)\cdot d}\f$ where \f$A_i(x) = a(x) \cdot 2^{-iK} \f$, that is, the polynomial corresponding to GLWEGadget limb \f$i\f$.

With a way to pack the required GLWEGadgets into GLWEs, the full protocol can be described below.

## Protocol description

The overall working of the protocol can be described as follows:

Preparation:

1. The client generates its secret key, and from it computes the "evaluation keys" for the server,
   specifically automorphism keys as will be required for query expansion and also a GGSW encryption
   of each component \f$s_i\f$ of the secret key \f$s\f$
2. The client sends the (public) evaluation keys to the server

Query access (repeat as many times as wanted):

1. Client: Use GLWEGadget packing to pack the GLWEGadget corresponding to the integer vector
           that is all 0 except in the position correponding to the index of the row to be selected
2. Client: Use GLWEGadget packing to pack the GLWEGadget corresponding to the integer vector
           corresponding to the binary representation of column number to be selected.
3. Client sends the 2 packed and encrypted queries to the server
4. Server: perform trace expansion into the row query
5. Server: Perform half-external products with the unpacked expanded values with each row of the database.
           Add each column's results together into a single GLWE
6. Server: query expand the column query and use external product to convert the resuling GLWEGadgets into
          GGSWs using the GGSW encryption(s) of \f$-s_i\f$
7. Server: Use a tree of CMux using as inputs the column sums and as selectors the expansion from the previous step
8. Server sends the resulting GLWE to the client
9. Client receives and decrypts result

The process is also pictured in the following figure:

\image html onionpircomp.svg "Representation of the OnionPIR protocol as in the example"
