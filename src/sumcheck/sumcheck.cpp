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
    // Prover's evaluation
    // Since zh(x) = x^l - 1, polynomialDivision is O(D)F
    vector<Fr> p = q; // Remainder
    p[0] -= s / l;

    vector<Fr> f = polynomialDivision(p, l);

    // Check if p is correct remainder i.e. degree < n-1 
    while (!p.empty() && p.back().isZero()) {
        p.pop_back();
    }

    if ((!p.empty() && p[0] != 0) || p.size() > l) throw runtime_error("Wrong remainder!");
    if (!p.empty()) p.erase(p.begin());
    else p.push_back(0);
    
    KZG::Commitment comm_f = commit(pk, f); // O(D)G
    KZG::Commitment comm_q = commit(pk, q); // O(D)G
    KZG::Commitment comm_p = commit(pk, p); // O(D)G

    // Prover sends comm_f and comm_q to Verifier
    // Verifier sends random challenge r to Prover
    Fr r;
    r.setByCSPRNG();

    // Prover creates witnesses to f(r) and q(r)
    KZG::Witness witness_f = createWitness(pk, f, r); // Witness to fr --> O(D)G
    KZG::Witness witness_q = createWitness(pk, q, r); // Witness to qr --> O(D)G
    KZG::Witness witness_p = createWitness(pk, p, r); // Witness to qr --> O(D)G

    // Prover sends to Verifier: witness_f, witness_q, witness_p

    // Public knowledge --> Both Prover and Verifier can evaluate Zh(r) in O(1)F
    Fr zr;
    Fr::pow(zr, r, l);
    zr -= 1;

    // Verifier's evaluation
    // V checks if the commitment and witness open to f(r) and q(r) --> O(1)G
    // V also checks that the evaluated qr = fr * zr --> O(1)G
    return verifyEval(pk, comm_f, r, witness_f) && verifyEval(pk, comm_q, r, witness_q) 
        && verifyEval(pk, comm_p, r, witness_p) && witness_q.qi == witness_f.qi * zr + s / l + r * witness_p.qi;
}
