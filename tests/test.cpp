#include "ntt.h"
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>

using namespace std;
using namespace mcl;
using namespace bn;

// Initialize the curve
void init_curve() {
    static bool initialized = false;
    if (!initialized) {
        mcl::bn::initPairing(mcl::BN_SNARK1);
        initialized = true;
    }
}

// Helper function to print a vector of Fr elements
void print_vector(const vector<Fr> &v, const string &name) {
    cout << name << ": [";
    for (size_t i = 0; i < v.size(); ++i) {
        cout << v[i];
        if (i < v.size() - 1) cout << ", ";
    }
    cout << "]" << endl;
}

// Test basic NTT functionality
void test_basic_ntt() {
    cout << "\n=== Testing Basic NTT ===" << endl;
    
    const int N = 8;
    vector<Fr> A(N);
    
    // Initialize with simple values
    for (int i = 0; i < N; ++i) {
        A[i] = i + 1;
    }
    
    print_vector(A, "Original");
    
    // Find primitive root
    Fr omega = findPrimitiveRoot(N);
    cout << "Primitive " << N << "-th root of unity: " << omega << endl;
    
    // Perform NTT
    vector<Fr> A_ntt = A;
    ntt_transform(A_ntt, omega);
    print_vector(A_ntt, "After NTT");
    
    // Verify that we can recover original with INTT
    // For now, we'll implement a simple inverse manually
    Fr omega_inv;
    Fr::inv(omega_inv, omega);
    vector<Fr> A_recovered = A_ntt;
    ntt_transform(A_recovered, omega_inv);
    
    // Scale by 1/N
    Fr n_inv;
    Fr::inv(n_inv, Fr(N));
    for (auto &x : A_recovered) {
        x *= n_inv;
    }
    
    print_vector(A_recovered, "After INTT");
    
    // Check if recovery is correct
    bool success = true;
    for (int i = 0; i < N; ++i) {
        if (A[i] != A_recovered[i]) {
            success = false;
            break;
        }
    }
    
    cout << "Recovery test: " << (success ? "PASSED" : "FAILED") << endl;
}

// Test with different sizes
void test_different_sizes() {
    cout << "\n=== Testing Different Sizes ===" << endl;
    
    vector<int> sizes = {4, 8, 16, 32};
    
    for (int N : sizes) {
        cout << "\nTesting size " << N << ":" << endl;
        
        vector<Fr> A(N);
        for (int i = 0; i < N; ++i) {
            A[i] = (i * i + 1) % 100; // Some pattern
        }
        
        try {
            Fr omega = findPrimitiveRoot(N);
            cout << "Found primitive root for N=" << N << endl;
            
            vector<Fr> A_copy = A;
            ntt_transform(A_copy, omega);
            cout << "NTT completed successfully" << endl;
            
        } catch (const exception &e) {
            cout << "Error: " << e.what() << endl;
        }
    }
}

// Test polynomial interpolation
void test_polynomial_interpolation() {
    cout << "\n=== Testing Polynomial Interpolation ===" << endl;
    
    const int N = 8;
    Fr omega = findPrimitiveRoot(N);
    
    // Original polynomial: p(x) = 1 + 2x + 3x^2
    vector<Fr> original_poly = {1, 2, 3, 0, 0, 0, 0, 0};
    print_vector(original_poly, "Original polynomial");
    
    // Evaluate polynomial at roots of unity (this is what NTT does)
    vector<Fr> evaluations = original_poly;
    ntt_transform(evaluations, omega);
    print_vector(evaluations, "Evaluations at roots of unity");
    
    // Now interpolate back using polynomial_interpolation
    // (This should give us back the original polynomial)
    vector<Fr> interpolated = evaluations;
    
    // Manual inverse NTT for interpolation
    Fr omega_inv;
    Fr::inv(omega_inv, omega);
    ntt_transform(interpolated, omega_inv);
    
    // Scale by 1/N
    Fr n_inv;
    Fr::inv(n_inv, Fr(N));
    for (auto &x : interpolated) {
        x *= n_inv;
    }
    
    print_vector(interpolated, "Interpolated polynomial");
    
    // Verify interpolation correctness
    bool success = true;
    for (int i = 0; i < N; ++i) {
        if (original_poly[i] != interpolated[i]) {
            success = false;
            cout << "Mismatch at coefficient " << i << ": " 
                 << original_poly[i] << " != " << interpolated[i] << endl;
            break;
        }
    }
    
    cout << "Interpolation test: " << (success ? "PASSED" : "FAILED") << endl;
}

// Test polynomial multiplication using NTT
void test_polynomial_multiply() {
    cout << "\n=== Testing Polynomial Multiplication ===" << endl;
    
    // Test case 1: Simple multiplication
    cout << "\nTest case 1: (1 + 2x) * (3 + x)" << endl;
    vector<Fr> p1 = {1, 2};        // 1 + 2x
    vector<Fr> p2 = {3, 1};        // 3 + x
    
    print_vector(p1, "Polynomial 1");
    print_vector(p2, "Polynomial 2");
    
    // Manual expected result: (1 + 2x)(3 + x) = 3 + x + 6x + 2x^2 = 3 + 7x + 2x^2
    vector<Fr> expected1 = {3, 7, 2};
    
    // Use NTT for multiplication
    int result_size = p1.size() + p2.size() - 1;
    int padded_size = 1;
    while (padded_size < result_size) padded_size <<= 1;
    
    // Pad polynomials
    p1.resize(padded_size, Fr(0));
    p2.resize(padded_size, Fr(0));
    
    Fr omega = findPrimitiveRoot(padded_size);
    
    // Transform
    vector<Fr> P1_ntt = p1, P2_ntt = p2;
    ntt_transform(P1_ntt, omega);
    ntt_transform(P2_ntt, omega);
    
    // Point-wise multiplication
    vector<Fr> result_ntt(padded_size);
    for (int i = 0; i < padded_size; ++i) {
        result_ntt[i] = P1_ntt[i] * P2_ntt[i];
    }
    
    // Inverse transform
    Fr omega_inv;
    Fr::inv(omega_inv, omega);
    ntt_transform(result_ntt, omega_inv);
    
    // Scale by 1/N
    Fr n_inv;
    Fr::inv(n_inv, Fr(padded_size));
    for (auto &x : result_ntt) {
        x *= n_inv;
    }
    
    // Trim to actual result size
    result_ntt.resize(result_size);
    print_vector(result_ntt, "NTT multiplication result");
    print_vector(expected1, "Expected result");
    
    // Verify
    bool success1 = (result_ntt.size() == expected1.size());
    if (success1) {
        for (size_t i = 0; i < expected1.size(); ++i) {
            if (result_ntt[i] != expected1[i]) {
                success1 = false;
                break;
            }
        }
    }
    cout << "Test case 1: " << (success1 ? "PASSED" : "FAILED") << endl;
    
    // Test case 2: Larger polynomials
    cout << "\nTest case 2: (1 + x + x^2) * (1 - x)" << endl;
    vector<Fr> p3 = {1, 1, 1};     // 1 + x + x^2
    vector<Fr> p4 = {1, -1};       // 1 - x
    
    print_vector(p3, "Polynomial 3");
    print_vector(p4, "Polynomial 4");
    
    // Manual expected: (1 + x + x^2)(1 - x) = 1 - x + x - x^2 + x^2 - x^3 = 1 - x^3
    vector<Fr> expected2 = {1, 0, 0, -1};
    
    // Prepare for NTT multiplication
    int result_size2 = p3.size() + p4.size() - 1;
    int padded_size2 = 1;
    while (padded_size2 < result_size2) padded_size2 <<= 1;
    
    p3.resize(padded_size2, Fr(0));
    p4.resize(padded_size2, Fr(0));
    
    Fr omega2 = findPrimitiveRoot(padded_size2);
    
    // Transform
    vector<Fr> P3_ntt = p3, P4_ntt = p4;
    ntt_transform(P3_ntt, omega2);
    ntt_transform(P4_ntt, omega2);
    
    // Point-wise multiplication
    vector<Fr> result2_ntt(padded_size2);
    for (int i = 0; i < padded_size2; ++i) {
        result2_ntt[i] = P3_ntt[i] * P4_ntt[i];
    }
    
    // Inverse transform
    Fr omega2_inv;
    Fr::inv(omega2_inv, omega2);
    ntt_transform(result2_ntt, omega2_inv);
    
    // Scale by 1/N
    Fr n2_inv;
    Fr::inv(n2_inv, Fr(padded_size2));
    for (auto &x : result2_ntt) {
        x *= n2_inv;
    }
    
    // Trim to actual result size
    result2_ntt.resize(result_size2);
    print_vector(result2_ntt, "NTT multiplication result");
    print_vector(expected2, "Expected result");
    
    // Verify
    bool success2 = (result2_ntt.size() == expected2.size());
    if (success2) {
        for (size_t i = 0; i < expected2.size(); ++i) {
            if (result2_ntt[i] != expected2[i]) {
                success2 = false;
                break;
            }
        }
    }
    cout << "Test case 2: " << (success2 ? "PASSED" : "FAILED") << endl;
}

// Test polynomial multiplication (conceptual) - kept for backward compatibility
void test_polynomial_concept() {
    cout << "\n=== Testing Polynomial Multiplication Concept ===" << endl;
    
    const int N = 8;
    
    // Two polynomials: p1(x) = 1 + 2x + 3x^2, p2(x) = 2 + x
    vector<Fr> p1 = {1, 2, 3, 0, 0, 0, 0, 0}; // Pad to size N
    vector<Fr> p2 = {2, 1, 0, 0, 0, 0, 0, 0}; // Pad to size N
    
    print_vector(p1, "Polynomial 1");
    print_vector(p2, "Polynomial 2");
    
    // Expected result: (1 + 2x + 3x^2)(2 + x) = 2 + x + 4x + 2x^2 + 6x^2 + 3x^3
    //                                           = 2 + 5x + 8x^2 + 3x^3
    
    Fr omega = findPrimitiveRoot(N);
    
    // Transform both polynomials
    vector<Fr> P1_ntt = p1, P2_ntt = p2;
    ntt_transform(P1_ntt, omega);
    ntt_transform(P2_ntt, omega);
    
    // Point-wise multiplication
    vector<Fr> result_ntt(N);
    for (int i = 0; i < N; ++i) {
        result_ntt[i] = P1_ntt[i] * P2_ntt[i];
    }
    
    // Inverse transform
    Fr omega_inv;
    Fr::inv(omega_inv, omega);
    ntt_transform(result_ntt, omega_inv);
    
    // Scale by 1/N
    Fr n_inv;
    Fr::inv(n_inv, Fr(N));
    for (auto &x : result_ntt) {
        x *= n_inv;
    }
    
    print_vector(result_ntt, "Product polynomial");
    
    // Manual verification
    cout << "Expected: [2, 5, 8, 3, 0, 0, 0, 0]" << endl;
}

// Test bit reversal function
void test_bit_reversal() {
    cout << "\n=== Testing Bit Reversal ===" << endl;
    
    int N = 8;
    int logN = 3;
    
    cout << "Bit reversal for N=" << N << " (logN=" << logN << "):" << endl;
    for (int i = 0; i < N; ++i) {
        int reversed = bitReverse(i, logN);
        cout << i << " -> " << reversed << endl;
    }
}

// Stress test with random data
void test_random_data() {
    cout << "\n=== Testing with Random Data ===" << endl;
    
    const int N = 16;
    vector<Fr> A(N);
    
    // Fill with random values
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(1, 1000);
    
    for (int i = 0; i < N; ++i) {
        A[i] = dis(gen);
    }
    
    cout << "Testing with " << N << " random values..." << endl;
    
    Fr omega = findPrimitiveRoot(N);
    
    // Store original
    vector<Fr> original = A;
    
    // Forward transform
    ntt_transform(A, omega);
    
    // Inverse transform
    Fr omega_inv;
    Fr::inv(omega_inv, omega);
    ntt_transform(A, omega_inv);
    
    // Scale by 1/N
    Fr n_inv;
    Fr::inv(n_inv, Fr(N));
    for (auto &x : A) {
        x *= n_inv;
    }
    
    // Check if we recovered the original
    bool success = true;
    for (int i = 0; i < N; ++i) {
        if (original[i] != A[i]) {
            success = false;
            cout << "Mismatch at position " << i << ": " 
                 << original[i] << " != " << A[i] << endl;
            break;
        }
    }
    
    cout << "Random data round-trip test: " << (success ? "PASSED" : "FAILED") << endl;
}

int main() {
    try {
        init_curve();
        cout << "MCL library initialized successfully" << endl;
        
        test_bit_reversal();
        test_basic_ntt();
        test_different_sizes();
        test_polynomial_interpolation();
        test_polynomial_multiply();
        test_polynomial_concept();
        test_random_data();
        
        cout << "\n=== All tests completed ===" << endl;
        
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
