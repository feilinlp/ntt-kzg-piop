#include "zerotest.h"
#include "../kzg/kzg.h"
#include "../ntt/ntt.h"
#include <mcl/bn.hpp>
#include <cassert>

using namespace std;
using namespace mcl;
using namespace bn;

// a = divided, b = divisor
vector<Fr> polynomialDivision(vector<Fr> a, vector<Fr> b) {
    // Remove leading zeros
    while (!a.empty() && a.back().isZero()) {
        a.pop_back();
    }
    while (!b.empty() && b.back().isZero()) {
        b.pop_back();
    }
    
    if (b.empty()) {
        throw runtime_error("Division by zero polynomial");
    }
    
    if (a.size() < b.size()) {
        return vector<Fr>(1, 0); // Result is zero
    }
    
    vector<Fr> quotient(a.size() - b.size() + 1, 0);
    vector<Fr> remainder = a;
    
    for (int i = quotient.size() - 1; i >= 0; i--) {
        if (remainder.size() >= b.size()) {
            // Calculate the coefficient for this term
            Fr coeff = remainder.back() / b.back();
            quotient[i] = coeff;
            
            // Subtract b * coeff * x^i from remainder
            for (int j = b.size() - 1; j >= 0; j--) {
                if (b[j] == 0) continue; // Reduces to O(1) if b = zh
                int pos = remainder.size() - b.size() + j;
                if (pos >= 0 && pos < remainder.size()) {
                    remainder[pos] = remainder[pos] - coeff * b[j];
                }
            }
            
            // Remove leading zero
            while (!remainder.empty() && remainder.back().isZero()) {
                remainder.pop_back();
            }
        }
    }
    
    return quotient;
}

// Proof size is O(1) as there is constant number of communication.
bool zeroTest(KZG::PublicKey pk, vector<Fr> q, Fr w, size_t l) {
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
    vector<Fr> f = polynomialDivision(q, zh);
    
    KZG::Commitment comm_f = commit(pk, f); // O(D)G
    KZG::Commitment comm_q = commit(pk, q); // O(D)G

    // Prover sends comm_f and comm_q to Verifier
    // Verifier sends random challenge r to Prover
    Fr r;
    r.setByCSPRNG();

    // Prover evaluates the values to f(r) and q(r)
    Fr qr = evaluatePolynomial(q, r);
    Fr fr = evaluatePolynomial(f, r);

    // Prover creates witnesses to f(r) and q(r)
    KZG::Witness witness_f = createWitness(pk, f, r); // Witness to fr --> O(D)G
    KZG::Witness witness_q = createWitness(pk, q, r); // Witness to qr --> O(D)G

    // Prover sends to Verifier: witness_f, witness_q, fr and qr

    // Public knowledge --> Both Prover and Verifier can evaluate Zh(r) in O(1)F
    Fr zr;
    Fr::pow(zr, r, l);
    zr -= 1;

    // Verifier's evaluation
    // V checks if the commitment and witness open to f(r) and q(r) --> O(1)G
    // V also checks that the evaluated qr = fr * zr --> O(1)G
    return verifyEval(pk, comm_f, r, fr, witness_f) && verifyEval(pk, comm_q, r, qr, witness_q) && qr == fr * zr;
}
