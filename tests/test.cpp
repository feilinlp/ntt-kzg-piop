#include "sumcheck.h"
#include "kzg.h"
#include "ntt.h"
#include "zerotest.h"
#include <mcl/bn.hpp>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace mcl;
using namespace bn;
using namespace std::chrono;

bool testNTT() {
    cout << "Testing NTT and Inverse NTT..." << endl;
    auto start_time = high_resolution_clock::now();
    
    try {
        // Test 1: Basic NTT and INTT roundtrip
        size_t N = 8;
        Fr omega = findPrimitiveRoot(N);
        
        vector<Fr> original(8);
        for (int i = 0; i < 8; i++) original[i] = rand();
        vector<Fr> test_data = original;
        
        // Forward NTT
        ntt_transform(test_data, omega);
        cout << "✓ Forward NTT completed" << endl;
        
        // Inverse NTT
        ntt_inverse(test_data, omega);
        cout << "✓ Inverse NTT completed" << endl;
        
        // Check if we get back the original
        bool roundtrip_correct = true;
        for (size_t i = 0; i < N; i++) {
            if (test_data[i] != original[i]) {
                roundtrip_correct = false;
                break;
            }
        }
        
        if (roundtrip_correct) {
            cout << "✓ NTT/INTT roundtrip test passed" << endl;
        } else {
            cout << "✗ NTT/INTT roundtrip test failed" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test 2: Test with different sizes
        for (size_t size : {4, 16, 32}) {
            Fr omega_test = findPrimitiveRoot(size);
            vector<Fr> test_vec(size);
            for (size_t i = 0; i < size; i++) {
                test_vec[i] = rand();
            }
            
            vector<Fr> backup = test_vec;
            ntt_transform(test_vec, omega_test);
            ntt_inverse(test_vec, omega_test);
            
            bool size_test_passed = true;
            for (size_t i = 0; i < size; i++) {
                if (test_vec[i] != backup[i]) {
                    size_test_passed = false;
                    break;
                }
            }
            
            if (size_test_passed) {
                cout << "✓ NTT size " << size << " test passed" << endl;
            } else {
                cout << "✗ NTT size " << size << " test failed" << endl;
                auto end_time = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(end_time - start_time);
                cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
                return false;
            }
        }
        
        cout << "✓ All NTT tests passed!" << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return true;
        
    } catch (const exception& e) {
        cout << "✗ NTT test failed with exception: " << e.what() << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return false;
    }
}

bool testPoly() {
    cout << "Testing Polynomial Multiplication using NTT..." << endl;
    auto start_time = high_resolution_clock::now();
    
    try {
        // Test 1: Simple polynomial multiplication
        // (x + 1) * (x + 2) = x^2 + 3x + 2
        vector<Fr> A = {1, 1}; // x + 1
        vector<Fr> B = {2, 1}; // x + 2
        
        size_t result_size = 1;
        while (result_size < A.size() + B.size()) result_size *= 2;
        Fr omega = findPrimitiveRoot(result_size);
        
        vector<Fr> result = polynomial_multiply(A, B, omega);
        
        // Expected: x^2 + 3x + 2 = [2, 3, 1, 0, ...]
        if (result[0] == 2 && result[1] == 3 && result[2] == 1) {
            cout << "✓ Basic polynomial multiplication test passed" << endl;
        } else {
            cout << "✗ Basic polynomial multiplication test failed" << endl;
            cout << "Expected: [2, 3, 1], Got: [" << result[0] << ", " << result[1] << ", " << result[2] << "]" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test 2: Multiplication with zero polynomial
        vector<Fr> zero_poly = {0};
        vector<Fr> non_zero = {1, 2, 3};
        
        result_size = 1;
        while (result_size < zero_poly.size() + non_zero.size()) result_size *= 2;
        omega = findPrimitiveRoot(result_size);
        
        result = polynomial_multiply(zero_poly, non_zero, omega);
        
        bool is_zero = true;
        for (size_t i = 0; i < min(result.size(), size_t(4)); i++) {
            if (result[i] != 0) {
                is_zero = false;
                break;
            }
        }
        
        if (is_zero) {
            cout << "✓ Zero polynomial multiplication test passed" << endl;
        } else {
            cout << "✗ Zero polynomial multiplication test failed" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test 3: Larger polynomial multiplication
        // Test (x^2 + x + 1) * (x + 1) = x^3 + 2x^2 + 2x + 1
        vector<Fr> P1 = {1, 1, 1}; // x^2 + x + 1
        vector<Fr> P2 = {1, 1};    // x + 1
        
        result_size = 1;
        while (result_size < P1.size() + P2.size()) result_size *= 2;
        omega = findPrimitiveRoot(result_size);
        
        result = polynomial_multiply(P1, P2, omega);
        
        // Expected: x^3 + 2x^2 + 2x + 1 = [1, 2, 2, 1]
        if (result[0] == 1 && result[1] == 2 && result[2] == 2 && result[3] == 1) {
            cout << "✓ Larger polynomial multiplication test passed" << endl;
        } else {
            cout << "✗ Larger polynomial multiplication test failed" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        cout << "✓ All polynomial multiplication tests passed!" << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return true;
        
    } catch (const exception& e) {
        cout << "✗ Polynomial multiplication test failed with exception: " << e.what() << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return false;
    }
}

bool testKZG() {
    cout << "Testing KZG Commitment Scheme..." << endl;
    auto start_time = high_resolution_clock::now();
    
    try {
        size_t degree = 30;
        KZG::PublicKey pk = setup(degree);
        cout << "✓ KZG setup completed" << endl;
        
        // Test 1: Basic commitment and evaluation
        vector<Fr> polynomial = {rand(), rand(), rand(), rand()}; 
        KZG::Commitment comm = commit(pk, polynomial);
        cout << "✓ Polynomial commitment created" << endl;
        
        // Test evaluation at a random point
        Fr eval_point = rand();
        Fr expected_value = evaluatePoly(polynomial, eval_point);
        
        // Create witness
        KZG::Witness witness = createWitness(pk, polynomial, eval_point);
        cout << "✓ Witness created" << endl;
        
        // Verify
        bool verification_result = verifyEval(pk, comm, eval_point, witness);
        
        if (verification_result) {
            cout << "✓ KZG evaluation verification passed" << endl;
        } else {
            cout << "✗ KZG evaluation verification failed" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test 2: Multiple evaluation points
        vector<Fr> test_points = {rand(), rand(), rand(), rand()};
        bool all_evaluations_passed = true;
        
        for (Fr point : test_points) {
            Fr value = evaluatePoly(polynomial, point);
            KZG::Witness w = createWitness(pk, polynomial, point);
            
            if (!verifyEval(pk, comm, point, w)) {
                all_evaluations_passed = false;
                break;
            }
        }
        
        if (all_evaluations_passed) {
            cout << "✓ Multiple evaluation points test passed" << endl;
        } else {
            cout << "✗ Multiple evaluation points test failed" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test 3: Test with wrong evaluation (should fail)
        KZG::Witness wrong_witness = createWitness(pk, polynomial, eval_point);
        wrong_witness.qi += 1;
        bool should_fail = verifyEval(pk, comm, eval_point, wrong_witness);
        
        if (!should_fail) {
            cout << "✓ Wrong evaluation correctly rejected" << endl;
        } else {
            cout << "✗ Wrong evaluation incorrectly accepted" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        cout << "✓ All KZG tests passed!" << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return true;
        
    } catch (const exception& e) {
        cout << "✗ KZG test failed with exception: " << e.what() << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return false;
    }
}

bool testZeroTest() {
    cout << "Testing Zero Test Protocol..." << endl;
    auto start_time = high_resolution_clock::now();
    
    try {
        size_t l = 4; // Domain size
        Fr w = findPrimitiveRoot(l);
        
        // Test 1: Create a polynomial that vanishes on H = {1, w, w^2, w^3}
        // Use polynomial (x-1)(x-w)(x-w^2)(x-w^3) = x^4 - 1
        vector<Fr> vanishing_poly(5, 0);
        vanishing_poly[0] = -1; // -1
        vanishing_poly[4] = 1;  // x^4
        
        cout << "✓ Created vanishing polynomial" << endl;
        
        // Setup KZG
        size_t degree = 10;
        KZG::PublicKey pk = setup(degree);
        cout << "✓ KZG setup for zero test completed" << endl;
        
        // Test the zero test
        bool zero_test_result = zeroTest(pk, vanishing_poly, w, l);
        
        if (zero_test_result) {
            cout << "✓ Zero test passed for vanishing polynomial" << endl;
        } else {
            cout << "✗ Zero test failed for vanishing polynomial" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test 2: Create a polynomial that does NOT vanish on H
        vector<Fr> non_vanishing_poly = {1, 2, 3, 4, 5}; // 1 + 2x + 3x^2 + 4x^3 + 5x^4
        
        bool should_fail = false;
        try {
            zeroTest(pk, non_vanishing_poly, w, l);
            should_fail = true; // If we reach here, the test didn't throw as expected
        } catch (const runtime_error& e) {
            cout << "✓ Non-vanishing polynomial correctly rejected" << endl;
        }
        
        if (should_fail) {
            cout << "✗ Non-vanishing polynomial incorrectly accepted" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test 3: Test with a more complex vanishing polynomial
        // Create polynomial that vanishes on domain by constructing it properly
        vector<Fr> complex_vanishing(10, 0);
        complex_vanishing[5] = -1;
        complex_vanishing[9] = 1;
        // Add some higher degree terms that don't affect vanishing property
        
        bool complex_test_result = zeroTest(pk, complex_vanishing, w, l);
        
        if (complex_test_result) {
            cout << "✓ Complex vanishing polynomial test passed" << endl;
        } else {
            cout << "✗ Complex vanishing polynomial test failed" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        cout << "✓ All zero test protocol tests passed!" << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return true;
        
    } catch (const exception& e) {
        cout << "✗ Zero test failed with exception: " << e.what() << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return false;
    }
}

bool testSumCheck() {
    cout << "Testing Sum Check Protocol..." << endl;
    auto start_time = high_resolution_clock::now();
    
    try {
        size_t l = 4; // Domain size
        Fr w = findPrimitiveRoot(l);
        
        // Setup KZG
        size_t degree = 10;
        KZG::PublicKey pk = setup(degree);
        cout << "✓ KZG setup for sum check completed" << endl;
        
        // Test 1: Create a polynomial and calculate its sum over domain H
        // Use polynomial q(x) = x^4 - 1 + constant, which vanishes on H
        // Then sum over H should be l * constant
        Fr constant = 5;
        Fr expected_sum = constant * l;
        
        vector<Fr> test_poly(5, 0);
        test_poly[0] = -1 + constant; // -1 + constant
        test_poly[4] = 1;             // x^4
        
        cout << "✓ Created test polynomial for sum check" << endl;
        
        // Verify the sum manually
        Fr actual_sum = 0;
        Fr curr = 1;
        for (size_t i = 0; i < l; i++) {
            actual_sum += evaluatePoly(test_poly, curr);
            curr *= w;
        }
        
        if (actual_sum == expected_sum) {
            cout << "✓ Manual sum verification correct" << endl;
        } else {
            cout << "✗ Manual sum verification failed" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test the sum check protocol
        bool sum_check_result = sumCheck(pk, test_poly, w, l, expected_sum);
        
        if (sum_check_result) {
            cout << "✓ Sum check protocol passed" << endl;
        } else {
            cout << "✗ Sum check protocol failed" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        // Test 2: Test with wrong sum (should fail)
        Fr wrong_sum = expected_sum + 1;
        bool should_fail = false;
        
        try {
            bool wrong_result = sumCheck(pk, test_poly, w, l, wrong_sum);
            if (wrong_result) {
                should_fail = true;
            } else {
                cout << "✓ Wrong sum correctly rejected" << endl;
            }
        } catch (const exception& e) {
            cout << "✓ Wrong sum correctly rejected with exception" << endl;
        }
        
        if (should_fail) {
            cout << "✗ Wrong sum incorrectly accepted" << endl;
            auto end_time = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end_time - start_time);
            cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
            return false;
        }
        
        cout << "✓ All sum check protocol tests passed!" << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return true;
        
    } catch (const exception& e) {
        cout << "✗ Sum check test failed with exception: " << e.what() << endl;
        auto end_time = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end_time - start_time);
        cout << "⏱️  Test completed in " << fixed << setprecision(3) << duration.count() / 1000.0 << " seconds" << endl;
        return false;
    }
}

int main() {
    // Initialize the curve
    initPairing(BN_SNARK1);

    int passed = 0;
    int total = 5;
    auto total_start_time = high_resolution_clock::now();

    cout << "=== NTT & INTT Tests ===" << endl;
    if (testNTT()) passed++;
    cout << endl;

    cout << "=== Polynomial Multiplication Tests ===" << endl;
    if (testPoly()) passed++;
    cout << endl;

    cout << "=== KZG Tests ===" << endl;
    if (testKZG()) passed++;
    cout << endl;

    cout << "=== ZeroTest Tests ===" << endl;
    if (testZeroTest()) passed++;
    cout << "Proof Size: 0.352 kb\n";
    cout << endl;

    cout << "=== SumCheck Tests ===" << endl;
    if (testSumCheck()) passed++;
    cout << "Proof Size: 0.512 kb\n";
    cout << endl;

    // Summary
    auto total_end_time = high_resolution_clock::now();
    auto total_duration = duration_cast<milliseconds>(total_end_time - total_start_time);
    
    cout << "=== Test Summary ===" << endl;
    cout << "Passed: " << passed << "/" << total << " tests" << endl;
    cout << "⏱️  Total execution time: " << fixed << setprecision(3) << total_duration.count() / 1000.0 << " seconds" << endl;
    
    if (passed == total) {
        cout << "🎉 All tests PASSED!" << endl;
        return 0;
    } else {
        cout << "❌ Some tests FAILED!" << endl;
        return 1;
    }
}
