#ifndef ZEROTEST_H
#define ZEROTEST_H

#include <mcl/bn.hpp>
#include "../kzg/kzg.h"
#include "../ntt/ntt.h"

using namespace mcl;
using namespace bn;
using namespace std;

vector<Fr> polynomialDivision(vector<Fr> &a, size_t n);

bool zeroTest(KZG::PublicKey pk, vector<Fr> q, Fr w, size_t l);

#endif // ZEROTEST_H
