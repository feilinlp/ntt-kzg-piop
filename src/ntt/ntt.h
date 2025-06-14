#ifndef NTT_H
#define NTT_H

#include <mcl/bn.hpp>
#include <vector>
#include <cassert>
#include <cmath>

using namespace mcl;
using namespace bn;
using namespace std;

/**
 * @brief Performs bit reversal for NTT
 * @param x The number to reverse bits for
 * @param logN Log base 2 of the array size
 * @return Bit-reversed number
 */
size_t bitReverse(size_t x, size_t logN);

/**
 * @brief Finds a primitive N-th root of unity in the finite field
 * @param N The order of the root (must be power of 2)
 * @return A primitive N-th root of unity
 */
Fr findPrimitiveRoot(size_t N);

/**
 * @brief Performs Number Theoretic Transform (NTT) on the input array
 * @param A Input/output vector of field elements (size must be power of 2)
 * @param omega Primitive N-th root of unity where N = A.size()
 * 
 * The function performs in-place NTT transformation.
 * After calling this function, A[i] will contain the i-th coefficient
 * of the NTT of the original polynomial.
 */
void ntt_transform(vector<Fr> &A, Fr omega);

/**
 * @brief Performs inverse Number Theoretic Transform (INTT)
 * @param A Input/output vector of field elements (size must be power of 2)
 * @param omega Primitive N-th root of unity where N = A.size()
 * 
 * This function undoes the NTT transformation by using the inverse
 * of omega and scaling by 1/N.
 */
void ntt_inverse(vector<Fr> &A, Fr omega);

/**
 * @brief Performs polynomial interpolation using inverse NTT
 * @param A Vector of polynomial evaluations at roots of unity
 * @param omega Primitive N-th root of unity where N = A.size()
 * @return Interpolated polynomial coefficients
 * 
 * Given evaluations of a polynomial at the N-th roots of unity
 * (omega^0, omega^1, ..., omega^(N-1)), this function recovers
 * the original polynomial coefficients using inverse NTT.
 * This is the inverse operation of polynomial evaluation at roots of unity.
 */
vector<Fr> polynomial_interpolation(vector<Fr> &A, Fr omega);

/**
 * @brief Multiplies two polynomials using NTT
 * @param A First polynomial coefficients
 * @param B Second polynomial coefficients
 * @return Product polynomial coefficients
 * 
 * The result size will be a.size() + b.size() - 1.
 * Input vectors will be zero-padded to the next power of 2 if needed.
 */
vector<Fr> polynomial_multiply(vector<Fr> &A, vector<Fr> &B, Fr omega);

#endif // NTT_H
