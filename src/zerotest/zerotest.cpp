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

    // Prover's evaluation
    // Vanishing Polynomial zh(x) = x^l - 1
    vector<Fr> zh(l+1, 0);
    zh[0] = -1;
    zh[l] = 1;

    // Since zh(x) = x^l - 1, polynomialDivision is O(D)F
    vector<Fr> f = polynomialDivision(q, zh);
    
    KZG::Commitment comm = commit(pk, f); // O(D)G

    // Prover sends comm to Verifier
    // Verifier sends random challenge r to Prover
    Fr r;
    r.setByCSPRNG();

    // Prover's evaluation
    Fr zr = evaluatePolynomial(zh, r); // O(D)F, same for fr and qr
    Fr fr = evaluatePolynomial(f, r); 
    Fr qr = evaluatePolynomial(q, r); 

    KZG::Witness witness = createWitness(pk, f, r); // Witness to fr --> O(D)G

    // Prover sends to Verifier: witness, qr, fr, zr

    // Verifier's evaluation
    // V checks if the commitment and witness open to f(r) --> O(1)G
    // V also checks that the evaluated qr equals to fr * zr --> O(1)F
    return verifyEval(pk, comm, witness) && qr == fr * zr;
}
