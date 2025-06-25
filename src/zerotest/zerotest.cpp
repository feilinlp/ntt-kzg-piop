#include "zerotest.h"
#include "../kzg/kzg.h"
#include "../ntt/ntt.h"
#include <mcl/bn.hpp>
#include <cassert>

using namespace std;
using namespace mcl;
using namespace bn;

vector<Fr> polynomialDivision(vector<Fr> &a, size_t n) {
    // Remove leading zeros
    while (!a.empty() && a.back().isZero()) {
        a.pop_back();
    }
    
    if (a.size() <= n) {
        return vector<Fr>(1, Fr(0));
    }
    
    vector<Fr> quotient(a.size() - n, Fr(0));
    
    // Copy high-degree coefficients to quotient and add them to low-degree terms
    for (int i = a.size() - 1; i >= n; i--) {
        quotient[i - n] = a[i];
        a[i - n] += a[i]; 
        a[i] = 0;
    }
    
    return quotient;
}

// Proof size is O(1) as there is constant number of communication.
bool zeroTest(KZG::PublicKey pk, vector<Fr> q, Fr w, size_t l) {
    // Check if q vanishes on H
    Fr curr = 1;
    for (size_t i = 0; i < l; i++) {
        if (evaluatePoly(q, curr) != 0) throw runtime_error("Polynomial does not vanish on H."); 
        curr *= w;
    }
    
    // Prover's evaluation
    // Since zh(x) = x^l - 1, polynomialDivision is O(D)F
    vector <Fr> remainder = q;
    vector<Fr> f = polynomialDivision(remainder, l);
    
    KZG::Commitment comm_f = commit(pk, f); // O(D)G
    KZG::Commitment comm_q = commit(pk, q); // O(D)G

    // Prover sends comm_f and comm_q to Verifier
    // Verifier sends random challenge r to Prover
    Fr r;
    r.setByCSPRNG();

    // Prover creates witnesses to f(r) and q(r)
    KZG::Witness witness_f = createWitness(pk, f, r); // Witness to fr --> O(D)G
    KZG::Witness witness_q = createWitness(pk, q, r); // Witness to qr --> O(D)G

    // Prover sends to Verifier: witness_f, witness_q

    // Public knowledge --> Both Prover and Verifier can evaluate Zh(r) in O(1)F
    Fr zr;
    Fr::pow(zr, r, l);
    zr -= 1;

    // Verifier's evaluation
    // V checks if the commitment and witness open to f(r) and q(r) --> O(1)G
    // V also checks that the evaluated qr = fr * zr --> O(1)G
    return verifyEval(pk, comm_f, r, witness_f) && verifyEval(pk, comm_q, r, witness_q) && witness_q.qi == witness_f.qi * zr;
}
