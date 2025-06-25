#ifndef KZG_H
#define KZG_H

#include <mcl/bn.hpp>

using namespace mcl;
using namespace bn;
using namespace std;

class KZG {
public:
    struct PublicKey {
        vector<G1> g1; 
        vector<G2> g2;
        size_t t; 
    };

    struct Commitment {
        G1 c;
    };

    struct Witness {
        Fr i;
        vector<Fr> q; // Polynomial (not needed)
        G1 w; // Witness
        Fr qi; // Evaluated value
    };
};

KZG::PublicKey setup(size_t t);

KZG::Commitment commit(KZG::PublicKey pk, vector<Fr> q);

Fr evaluatePoly(vector<Fr> q, Fr i);

vector<Fr> divideByLinear(vector<Fr> q, Fr i);

KZG::Witness createWitness(KZG::PublicKey pk, vector<Fr> q, Fr i);

bool verifyEval(KZG::PublicKey pk, KZG::Commitment comm, Fr i, KZG::Witness witness);

#endif // KZG_H
