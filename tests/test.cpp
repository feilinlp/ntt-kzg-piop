#include "ntt.h"
#include "kzg.h"
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
        initPairing(BN_SNARK1);
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

// Test KZG setup
void test_kzg_setup() {
    cout << "\n=== Testing KZG Setup ===" << endl;
    
    size_t t = 10;
    auto pk = setup(t);
    
    cout << "Setup completed for degree " << t << endl;
    cout << "G1 points generated: " << pk.g1.size() << endl;
    cout << "G2 points generated: " << pk.g2.size() << endl;
    cout << "Parameter t: " << pk.t << endl;
    
    // Verify we have the correct number of points
    bool size_correct = (pk.g1.size() == t + 1 && pk.g2.size() == t + 1);
    cout << "Size check: " << (size_correct ? "PASSED" : "FAILED") << endl;
    
    // Verify points are not zero (they should be different from identity)
    bool non_zero = true;
    G1 zero_g1;
    G2 zero_g2;
    
    for (size_t i = 0; i <= t; ++i) {
        if (pk.g1[i] == zero_g1 || pk.g2[i] == zero_g2) {
            non_zero = false;
            cout << "Found zero point at index " << i << endl;
            break;
        }
    }
    cout << "Non-zero points check: " << (non_zero ? "PASSED" : "FAILED") << endl;
}

// Test polynomial commitment
void test_kzg_commit() {
    cout << "\n=== Testing KZG Commitment ===" << endl;
    
    size_t t = 5;
    auto pk = setup(t);
    
    // Test polynomial: 1 + 2x + 3x^2
    vector<Fr> poly = {1, 2, 3};
    
    cout << "Polynomial coefficients: [";
    for (size_t i = 0; i < poly.size(); ++i) {
        cout << poly[i];
        if (i < poly.size() - 1) cout << ", ";
    }
    cout << "]" << endl;
    
    auto commitment = commit(pk, poly);
    
    // Verify commitment is not zero
    G1 zero_g1;
    bool non_zero_commit = !(commitment.c == zero_g1);
    cout << "Non-zero commitment check: " << (non_zero_commit ? "PASSED" : "FAILED") << endl;
    
    // Test empty polynomial commitment
    vector<Fr> empty_poly;
    auto empty_commit = commit(pk, empty_poly);
    bool empty_is_zero = (empty_commit.c == zero_g1);
    cout << "Empty polynomial commitment is zero: " << (empty_is_zero ? "PASSED" : "FAILED") << endl;
    
    // Test constant polynomial
    vector<Fr> const_poly = {5};
    auto const_commit = commit(pk, const_poly);
    bool const_non_zero = !(const_commit.c == zero_g1);
    cout << "Constant polynomial commitment non-zero: " << (const_non_zero ? "PASSED" : "FAILED") << endl;
}

// Test polynomial evaluation
void test_polynomial_evaluation() {
    cout << "\n=== Testing Polynomial Evaluation ===" << endl;
    
    // Test polynomial: 1 + 2x + 3x^2
    vector<Fr> poly = {1, 2, 3};
    
    // Test evaluation at x = 0: should be 1
    Fr x0 = 0;
    Fr result0 = evaluatePolynomial(poly, x0);
    bool test0 = (result0 == Fr(1));
    cout << "p(0) = " << result0 << ", expected 1: " << (test0 ? "PASSED" : "FAILED") << endl;
    
    // Test evaluation at x = 1: should be 1 + 2 + 3 = 6
    Fr x1 = 1;
    Fr result1 = evaluatePolynomial(poly, x1);
    bool test1 = (result1 == Fr(6));
    cout << "p(1) = " << result1 << ", expected 6: " << (test1 ? "PASSED" : "FAILED") << endl;
    
    // Test evaluation at x = 2: should be 1 + 4 + 12 = 17
    Fr x2 = 2;
    Fr result2 = evaluatePolynomial(poly, x2);
    bool test2 = (result2 == Fr(17));
    cout << "p(2) = " << result2 << ", expected 17: " << (test2 ? "PASSED" : "FAILED") << endl;
    
    // Test with single coefficient
    vector<Fr> single = {42};
    Fr single_result = evaluatePolynomial(single, Fr(100));
    bool single_test = (single_result == Fr(42));
    cout << "Constant polynomial test: " << (single_test ? "PASSED" : "FAILED") << endl;
}

// Test witness creation and verification
void test_kzg_witness() {
    cout << "\n=== Testing KZG Witness Creation and Verification ===" << endl;
    
    size_t t = 10;
    auto pk = setup(t);
    
    // Test polynomial: 1 + 2x + 3x^2 + x^3
    vector<Fr> poly = {1, 2, 3, 1};
    auto commitment = commit(pk, poly);
    
    cout << "Testing polynomial: 1 + 2x + 3x^2 + x^3" << endl;
    
    // Test witness creation and verification at different points
    vector<Fr> test_points = {0, 1, 2, 5, 10};
    
    for (const auto& point : test_points) {
        cout << "\nTesting at x = " << point << ":" << endl;
        
        // Create witness
        auto witness = createWitness(pk, poly, point);
        
        // Verify the witness
        bool verification = verifyEval(pk, commitment, witness);
        cout << "Verification result: " << (verification ? "PASSED" : "FAILED") << endl;
        
        // Double-check by computing expected value
        Fr expected_value = evaluatePolynomial(poly, point);
        Fr witness_value = evaluatePolynomial(witness.q, witness.i);
        cout << "Expected p(" << point << ") = " << expected_value << endl;
        cout << "Witness stores p(" << witness.i << ") = " << witness_value << endl;
    }
}

// Test edge cases for KZG
void test_kzg_edge_cases() {
    cout << "\n=== Testing KZG Edge Cases ===" << endl;
    
    size_t t = 3;
    auto pk = setup(t);
    
    // Test 1: Single coefficient polynomial
    cout << "\nTest 1: Constant polynomial" << endl;
    vector<Fr> const_poly = {7};
    auto const_commit = commit(pk, const_poly);
    auto const_witness = createWitness(pk, const_poly, Fr(5));
    bool const_verify = verifyEval(pk, const_commit, const_witness);
    cout << "Constant polynomial verification: " << (const_verify ? "PASSED" : "FAILED") << endl;
    
    // Test 2: Linear polynomial
    cout << "\nTest 2: Linear polynomial" << endl;
    vector<Fr> linear_poly = {3, 4}; // 3 + 4x
    auto linear_commit = commit(pk, linear_poly);
    auto linear_witness = createWitness(pk, linear_poly, Fr(2));
    bool linear_verify = verifyEval(pk, linear_commit, linear_witness);
    cout << "Linear polynomial verification: " << (linear_verify ? "PASSED" : "FAILED") << endl;
    
    // Test 3: Zero polynomial
    cout << "\nTest 3: Zero polynomial" << endl;
    vector<Fr> zero_poly = {0};
    auto zero_commit = commit(pk, zero_poly);
    auto zero_witness = createWitness(pk, zero_poly, Fr(1));
    bool zero_verify = verifyEval(pk, zero_commit, zero_witness);
    cout << "Zero polynomial verification: " << (zero_verify ? "PASSED" : "FAILED") << endl;
    
    // Test 4: Maximum degree polynomial
    cout << "\nTest 4: Maximum degree polynomial" << endl;
    vector<Fr> max_poly(t + 1);
    for (size_t i = 0; i <= t; ++i) {
        max_poly[i] = i + 1; // 1 + 2x + 3x^2 + ... + (t+1)x^t
    }
    auto max_commit = commit(pk, max_poly);
    auto max_witness = createWitness(pk, max_poly, Fr(1));
    bool max_verify = verifyEval(pk, max_commit, max_witness);
    cout << "Maximum degree polynomial verification: " << (max_verify ? "PASSED" : "FAILED") << endl;
}

// Test multiple witnesses for same polynomial
void test_kzg_multiple_witnesses() {
    cout << "\n=== Testing Multiple Witnesses for Same Polynomial ===" << endl;
    
    size_t t = 8;
    auto pk = setup(t);
    
    // Test polynomial: x^3 - 2x^2 + x - 1
    vector<Fr> poly = {-1, 1, -2, 1};
    auto commitment = commit(pk, poly);
    
    cout << "Testing polynomial: x^3 - 2x^2 + x - 1" << endl;
    
    // Create multiple witnesses
    vector<Fr> points = {0, 1, 2, 3, -1};
    vector<KZG::Witness> witnesses;
    
    for (const auto& point : points) {
        witnesses.push_back(createWitness(pk, poly, point));
    }
    
    // Verify all witnesses
    bool all_passed = true;
    for (size_t i = 0; i < witnesses.size(); ++i) {
        bool result = verifyEval(pk, commitment, witnesses[i]);
        cout << "Witness " << i << " (x=" << points[i] << "): " << (result ? "PASSED" : "FAILED") << endl;
        if (!result) all_passed = false;
    }
    
    cout << "All witnesses verification: " << (all_passed ? "PASSED" : "FAILED") << endl;
}

// Test with random polynomials
void test_kzg_random() {
    cout << "\n=== Testing KZG with Random Polynomials ===" << endl;
    
    size_t t = 15;
    auto pk = setup(t);
    
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(1, 100);
    
    // Generate random polynomial
    size_t degree = 5 + (dis(gen) % 6); // degree between 5 and 10
    vector<Fr> poly(degree + 1);
    
    cout << "Random polynomial of degree " << degree << ": [";
    for (size_t i = 0; i <= degree; ++i) {
        poly[i] = dis(gen);
        cout << poly[i];
        if (i < degree) cout << ", ";
    }
    cout << "]" << endl;
    
    auto commitment = commit(pk, poly);
    
    // Test at random points
    int num_tests = 5;
    bool all_passed = true;
    
    for (int i = 0; i < num_tests; ++i) {
        Fr test_point = dis(gen);
        auto witness = createWitness(pk, poly, test_point);
        bool result = verifyEval(pk, commitment, witness);
        
        cout << "Test " << (i+1) << " at x=" << test_point << ": " << (result ? "PASSED" : "FAILED") << endl;
        if (!result) all_passed = false;
    }
    
    cout << "Random polynomial tests: " << (all_passed ? "PASSED" : "FAILED") << endl;
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
        
        // KZG Tests
        test_kzg_setup();
        test_kzg_commit();
        test_polynomial_evaluation();
        test_kzg_witness();
        test_kzg_edge_cases();
        test_kzg_multiple_witnesses();
        test_kzg_random();
        
        cout << "\n=== All tests completed ===" << endl;
        
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
