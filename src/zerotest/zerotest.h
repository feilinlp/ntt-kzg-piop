#ifndef ZEROTEST_H
#define ZEROTEST_H

#include <mcl/bn.hpp>
#include <chrono>
#include <iomanip>
#include "../kzg/kzg.h"
#include "../ntt/ntt.h"

using namespace mcl;
using namespace bn;
using namespace std;
using namespace std::chrono;

vector<Fr> polynomialDivision(vector<Fr> &a, size_t n);

void startTime(high_resolution_clock::time_point &start_time);

void endTime(milliseconds &time, high_resolution_clock::time_point &start_time);

bool zeroTest(KZG::PublicKey pk, vector<Fr> q, Fr w, size_t l);

#endif // ZEROTEST_H
