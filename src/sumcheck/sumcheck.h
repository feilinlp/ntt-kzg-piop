#ifndef SUMCHECK_H
#define SUMCHECK_H

#include <mcl/bn.hpp>
#include "../kzg/kzg.h"
#include "../ntt/ntt.h"

using namespace mcl;
using namespace bn;
using namespace std;

bool sumCheck(KZG::PublicKey pk, vector<Fr> q, Fr w, size_t l, Fr s);

#endif // SUMCHECK_H
