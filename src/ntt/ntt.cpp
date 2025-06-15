#include <mcl/bn.hpp>
#include <iostream>
#include <cmath>

using namespace mcl;
using namespace bn;
using namespace std;

size_t bitReverse(size_t x, size_t logN) {
    size_t res = 0;
    for (size_t i = 0; i < logN; ++i) {
        res <<= 1;
        res |= (x >> i) & 1;
    }
    return res;
}

Fr findPrimitiveRoot(size_t N) {
    Fr order = -1;

    Fr g, omega, tmp;
    for (size_t i = 2; ; ++i) {
        g = i;
        Fr::pow(omega, g, order / N);

        // Check: omega^N == 1
        Fr::pow(tmp, omega, N);
        if (tmp != 1) continue;

        // Check: omega^(N/2) == -1 
        Fr::pow(tmp, omega, N / 2);
        if (tmp != -1) continue;

        return omega;
    }
}

void ntt_transform(vector<Fr> &A, Fr omega) {
    size_t n = A.size();
    size_t logN = log2(n);
    assert((1 << logN) == n); // n must be power of 2

    for (size_t i = 0; i < n; ++i) {
        size_t j = bitReverse(i, logN);
        if (i < j) swap(A[i], A[j]);
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        Fr wlen;
        Fr::pow(wlen, omega, n / len);

        for (size_t i = 0; i < n; i += len) {
            Fr w = 1;
            for (size_t j = 0; j < len / 2; ++j) {
                Fr u = A[i + j];
                Fr v = A[i + j + len / 2] * w;
                A[i + j] = u + v;
                A[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

void ntt_inverse(vector<Fr> &A, Fr omega) {
    size_t n = A.size();

    Fr omega_inv;
    Fr::inv(omega_inv, omega);
    ntt_transform(A, omega_inv);

    Fr n_inv;
    Fr::inv(n_inv, n);
    for (auto &x : A) {
        x *= n_inv;
    }
}

vector<Fr> polynomial_interpolation(vector<Fr> &A, Fr omega) {
    ntt_inverse(A, omega);
    return A;
}

vector<Fr> polynomial_multiply(vector<Fr> &A, vector<Fr> &B, Fr omega) {
    size_t A_n = A.size();
    size_t B_n = B.size();

    size_t n = 1;
    while (n < A_n + B_n) n *= 2;

    A.resize(n, 0);
    B.resize(n, 0);

    ntt_transform(A, omega);
    ntt_transform(B, omega);

    vector<Fr> result(n);
    for (size_t i = 0; i < n; i++) {
        result[i] = A[i] * B[i];
    }

    ntt_inverse(result, omega);
    return result;
}
