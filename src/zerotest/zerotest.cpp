#include "zerotest.h"
#include "../kzg/kzg.h"
#include "../ntt/ntt.h"
#include <mcl/bn.hpp>
#include <cassert>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace mcl;
using namespace bn;
using namespace std::chrono;

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

void startTime(high_resolution_clock::time_point &start_time) {
    start_time = high_resolution_clock::now();
}

void endTime(milliseconds &time, high_resolution_clock::time_point &start_time) {
    auto end_time = high_resolution_clock::now();
    time += duration_cast<milliseconds>(end_time - start_time);
}

// Proof size is O(1) as there is constant number of communication.
bool zeroTest(KZG::PublicKey pk, vector<Fr> q, Fr w, size_t l) {    
    auto start_time = high_resolution_clock::now();
    
    milliseconds prover_time = duration_cast<milliseconds>(start_time - start_time);
    milliseconds verifier_time = duration_cast<milliseconds>(start_time - start_time);

    // Check if q vanishes on H
    Fr curr = 1;
    for (size_t i = 0; i < l; i++) {
        if (evaluatePoly(q, curr) != 0) throw runtime_error("Polynomial does not vanish on H."); 
        curr *= w;
    }
    
    // Prover's evaluation
    // Since zh(x) = x^l - 1, polynomialDivision is O(D)F
    startTime(start_time);
    vector <Fr> remainder = q;
    vector<Fr> f = polynomialDivision(remainder, l);
    
    KZG::Commitment comm_f = commit(pk, f); // O(D)G
    KZG::Commitment comm_q = commit(pk, q); // O(D)G
    endTime(prover_time, start_time);

    // Prover sends comm_f and comm_q to Verifier
    // Verifier sends random challenge r to Prover
    startTime(start_time);
    Fr r;
    r.setByCSPRNG();
    endTime(verifier_time, start_time);

    // Prover creates witnesses to f(r) and q(r)
    startTime(start_time);
    KZG::Witness witness_f = createWitness(pk, f, r); // Witness to fr --> O(D)G
    KZG::Witness witness_q = createWitness(pk, q, r); // Witness to qr --> O(D)G
    endTime(prover_time, start_time);

    // Prover sends to Verifier: witness_f, witness_q

    // Public knowledge --> Both Prover and Verifier can evaluate Zh(r) in O(1)F
    startTime(start_time);
    Fr zr;
    Fr::pow(zr, r, l);
    zr -= 1;

    // Verifier's evaluation
    // V checks if the commitment and witness open to f(r) and q(r) --> O(1)G
    // V also checks that the evaluated qr = fr * zr --> O(1)G
    bool succeed = verifyEval(pk, comm_f, r, witness_f) && verifyEval(pk, comm_q, r, witness_q) && witness_q.qi == witness_f.qi * zr;
    endTime(verifier_time, start_time);

    cout << "\nRunning ZeroTest...\n";
    cout << "Prover time: " << fixed << setprecision(3) << prover_time.count() / 1000.0 << " seconds" << endl;
    cout << "Verifier time: " << fixed << setprecision(3) << verifier_time.count() / 1000.0 << " seconds" << endl;

    return succeed;
}
