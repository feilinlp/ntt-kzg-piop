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

    // Vanishing Polynomial zh(x) = x^l - 1
    vector<Fr> zh(l+1, 0);
    zh[0] = -1;
    zh[l] = 1;
    
    // Prover's evaluation
    // Since zh(x) = x^l - 1, polynomialDivision is O(D)F
    vector<Fr> temp = q;
    temp[0] -= s / l;
    vector<Fr> f = polynomialDivision(temp, zh);
    
    KZG::Commitment comm_f = commit(pk, f); // O(D)G
    KZG::Commitment comm_q = commit(pk, q); // O(D)G

    // Prover sends comm_f and comm_q to Verifier
    // Verifier sends random challenge r to Prover
    Fr r;
    r.setByCSPRNG();

    // Prover evaluates the values to f(r) and q(r)
    KZG::Witness witness_f = createWitness(pk, f, r); // Witness to fr --> O(D)G
    KZG::Witness witness_q = createWitness(pk, q, r); // Witness to qr --> O(D)G

    // Prover sends to Verifier: witness_f and witness_q

    // Public knowledge --> Both Prover and Verifier can evaluate Zh(r) in O(1)F
    Fr zr;
    Fr::pow(zr, r, l);
    zr -= 1;

    // Verifier's evaluation
    // V checks if the commitment and witness open to f(r) and q(r) --> O(1)G
    // V also checks that the evaluated qr = fr * zr --> O(1)F
    return verifyEval(pk, comm_f, witness_f) && verifyEval(pk, comm_q, witness_q) && witness_q.w == witness_f.w * zr + s / l;
}
