#include "sumcheck.h"
#include "../zerotest/zerotest.h"
#include "../kzg/kzg.h"
#include "../ntt/ntt.h"
#include <mcl/bn.hpp>
#include <cassert>

using namespace std;
using namespace mcl;
using namespace bn;

void outputTiming(milliseconds prover_time, milliseconds verifier_time) {
    cout << "\nRunning SumCheck...\n";
    cout << "Prover time: " << fixed << setprecision(3) << prover_time.count() / 1000.0 << " seconds" << endl;
    cout << "Verifier time: " << fixed << setprecision(3) << verifier_time.count() / 1000.0 << " seconds" << endl;
}

// Proof size is O(1) as there is constant number of communication.
bool sumCheck(KZG::PublicKey pk, vector<Fr> q, Fr w, size_t l, Fr s) {  
    auto start_time = high_resolution_clock::now();

    milliseconds prover_time = duration_cast<milliseconds>(start_time - start_time);
    milliseconds verifier_time = duration_cast<milliseconds>(start_time - start_time);  

    // Prover's evaluation
    // Since zh(x) = x^l - 1, polynomialDivision is O(D)F
    startTime(start_time);
    vector<Fr> p = q; // Remainder
    p[0] -= s / l;

    vector<Fr> f = polynomialDivision(p, l);

    // Check if p is correct remainder i.e. degree < n-1 
    while (!p.empty() && p.back().isZero()) {
        p.pop_back();
    }

    if ((!p.empty() && p[0] != 0) || p.size() > l) {
        outputTiming(prover_time, verifier_time); 
        throw runtime_error("Wrong remainder!");
    }

    if (!p.empty()) p.erase(p.begin());
    else p.push_back(0);
    
    KZG::Commitment comm_f = commit(pk, f); // O(D)G
    KZG::Commitment comm_q = commit(pk, q); // O(D)G
    KZG::Commitment comm_p = commit(pk, p); // O(D)G
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
    KZG::Witness witness_p = createWitness(pk, p, r); // Witness to qr --> O(D)G
    endTime(prover_time, start_time);

    // Prover sends to Verifier: witness_f, witness_q, witness_p

    // Public knowledge --> Both Prover and Verifier can evaluate Zh(r) in O(1)F
    startTime(start_time);
    Fr zr;
    Fr::pow(zr, r, l);
    zr -= 1;
    
    // Verifier's evaluation
    // V checks if the commitment and witness open to f(r) and q(r) --> O(1)G
    // V also checks that the evaluated qr = fr * zr --> O(1)G
    bool succeed = verifyEval(pk, comm_f, r, witness_f) && verifyEval(pk, comm_q, r, witness_q) 
        && verifyEval(pk, comm_p, r, witness_p) && witness_q.qi == witness_f.qi * zr + s / l + r * witness_p.qi;
    endTime(verifier_time, start_time);
    
    outputTiming(prover_time, verifier_time); 

    return succeed;
}
