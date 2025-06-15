#include "sumcheck.h"
#include "../zerotest/zerotest.h"
#include "../kzg/kzg.h"
#include "../ntt/ntt.h"
#include <mcl/bn.hpp>
#include <cassert>

using namespace std;
using namespace mcl;
using namespace bn;

// Proof size is O(1) as there is constant number of communication.
bool sumCheck(KZG::PublicKey pk, vector<Fr> q, Fr w, size_t l, Fr s) {
    // Check if q vanishes on H
    Fr curr = 1;
    for (size_t i = 0; i < l; i++) {
        if (evaluatePolynomial(q, curr) != 0) throw runtime_error("Polynomial does not vanish on H."); 
        curr *= w;
    }

    // Prover's evaluation
    // Vanishing Polynomial zh(x) = x^l - 1
    vector<Fr> zh(l+1, 0);
    zh[0] = -1;
    zh[l] = 1;

    // Since zh(x) = x^l - 1, polynomialDivision is O(D)F
    vector<Fr> temp = q;
    temp[0] -= s / l;
    vector<Fr> f = polynomialDivision(temp, zh);
    
    KZG::Commitment comm = commit(pk, f); // O(D)G

    // Prover sends s and comm to Verifier
    // Verifier sends random challenge r to Prover
    Fr r;
    r.setByCSPRNG();

    // Prover's evaluation
    Fr zr = evaluatePolynomial(zh, r); // O(D)F, same for fr and qr
    Fr fr = evaluatePolynomial(f, r); 
    Fr qr = evaluatePolynomial(q, r); 

    // What happens after this???
    KZG::Witness witness = createWitness(pk, f, r); // Witness to fr --> O(D)G

    // Prover sends to Verifier: witness, qr, fr, zr

    // Verifier's evaluation
    // V checks if the commitment and witness open to f(r) --> O(1)G
    // V also checks that the evaluated qr = fr * zr + s / l --> O(1)F
    return verifyEval(pk, comm, witness) && qr == fr * zr + s / l;
}
