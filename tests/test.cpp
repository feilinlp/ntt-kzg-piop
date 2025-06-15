#include "sumcheck.h"
#include "kzg.h"
#include "ntt.h"
#include <mcl/bn.hpp>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

using namespace std;
using namespace mcl;
using namespace bn;

// Helper function to check if n is a power of 2
bool isPowerOfTwo(size_t n) {
    return n > 0 && (n & (n - 1)) == 0;
}

// Helper function to create a polynomial that vanishes on H and has a specific sum
vector<Fr> createValidSumCheckPolynomial(size_t l, Fr omega, Fr target_sum) {
    assert(isPowerOfTwo(l));
    
    // Create vanishing polynomial zh(x) = x^l - 1
    vector<Fr> zh(l + 1, 0);
    zh[0] = -1;  // -1
    zh[l] = 1;   // x^l
    
    // Create a random polynomial f(x) of degree < l
    vector<Fr> f(l);
    for (size_t i = 0; i < l; i++) {
        f[i].setByCSPRNG();
    }
    
    // Compute f * zh using polynomial multiplication
    vector<Fr> f_copy = f;
    vector<Fr> zh_copy = zh;
    
    // For polynomial multiplication, we need a larger domain
    size_t mult_size = 1;
    while (mult_size < f.size() + zh.size()) mult_size *= 2;
    
    Fr mult_omega = findPrimitiveRoot(mult_size);
    vector<Fr> fzh = polynomial_multiply(f_copy, zh_copy, mult_omega);
    
    // Resize to appropriate degree
    fzh.resize(l + f.size());
    
    // Add constant term to make sum equal target_sum
    // Since f*zh vanishes on H, we just need to add target_sum/l
    Fr constant_term = target_sum;
    Fr l_inv;
    Fr l_fr = l;
    Fr::inv(l_inv, l_fr);
    Fr::mul(constant_term, constant_term, l_inv);
    
    fzh[0] = fzh[0] + constant_term;
    
    return fzh;
}

// Test case 1: Basic valid sum check with power-of-2 domain
bool testBasicValidSumCheck() {
    cout << "Test 1: Basic Valid Sum Check (l=4)..." << endl;
    
    try {
        size_t l = 4;  // Must be power of 2 for NTT
        KZG::PublicKey pk = setup(3 * l);  // Larger setup for safety
        
        // Find primitive root using NTT function
        Fr omega = findPrimitiveRoot(l);
        
        // Target sum
        Fr target_sum = 42;
        
        // Create valid polynomial
        vector<Fr> q = createValidSumCheckPolynomial(l, omega, target_sum);
        
        bool result = sumCheck(pk, q, omega, l, target_sum);
        
        if (result) {
            cout << "✓ Test 1 PASSED" << endl;
            return true;
        } else {
            cout << "✗ Test 1 FAILED" << endl;
            return false;
        }
    } catch (const exception& e) {
        cout << "✗ Test 1 FAILED with exception: " << e.what() << endl;
        return false;
    }
}

// Test case 2: Invalid sum should fail
bool testInvalidSum() {
    cout << "Test 2: Invalid Sum Check..." << endl;
    
    try {
        size_t l = 4;
        KZG::PublicKey pk = setup(3 * l);
        
        Fr omega = findPrimitiveRoot(l);
        
        Fr correct_sum = 100;
        Fr wrong_sum = 200;
        
        // Create polynomial that sums to correct_sum
        vector<Fr> q = createValidSumCheckPolynomial(l, omega, correct_sum);
        
        // But claim it sums to wrong_sum
        bool result = sumCheck(pk, q, omega, l, wrong_sum);
        
        if (!result) {
            cout << "✓ Test 2 PASSED (correctly rejected invalid sum)" << endl;
            return true;
        } else {
            cout << "✗ Test 2 FAILED (should have rejected invalid sum)" << endl;
            return false;
        }
    } catch (const exception& e) {
        cout << "✗ Test 2 FAILED with exception: " << e.what() << endl;
        return false;
    }
}

// Test case 3: Non-vanishing polynomial should throw exception
bool testNonVanishingPolynomial() {
    cout << "Test 3: Non-Vanishing Polynomial..." << endl;
    
    try {
        size_t l = 4;
        KZG::PublicKey pk = setup(3 * l);
        
        Fr omega = findPrimitiveRoot(l);
        
        // Create a polynomial that does NOT vanish on H
        // Simple constant polynomial
        vector<Fr> q(5, 0);
        q[0] = 1;  // q(x) = 1, which doesn't vanish anywhere
        
        Fr s = l;  // This would be the sum if it were valid
        
        // This should throw an exception
        bool result = sumCheck(pk, q, omega, l, s);
        
        cout << "✗ Test 3 FAILED (should have thrown exception)" << endl;
        return false;
        
    } catch (const runtime_error& e) {
        cout << "✓ Test 3 PASSED (correctly caught non-vanishing polynomial)" << endl;
        return true;
    } catch (const exception& e) {
        cout << "✗ Test 3 FAILED with unexpected exception: " << e.what() << endl;
        return false;
    }
}

// Test case 4: Larger domain (l=8)
bool testLargerDomain() {
    cout << "Test 4: Larger Domain (l=8)..." << endl;
    
    try {
        size_t l = 8;
        KZG::PublicKey pk = setup(3 * l);
        
        Fr omega = findPrimitiveRoot(l);
        
        Fr target_sum = 1337;
        vector<Fr> q = createValidSumCheckPolynomial(l, omega, target_sum);
        
        bool result = sumCheck(pk, q, omega, l, target_sum);
        
        if (result) {
            cout << "✓ Test 4 PASSED" << endl;
            return true;
        } else {
            cout << "✗ Test 4 FAILED" << endl;
            return false;
        }
    } catch (const exception& e) {
        cout << "✗ Test 4 FAILED with exception: " << e.what() << endl;
        return false;
    }
}

// Test case 5: Zero sum
bool testZeroSum() {
    cout << "Test 5: Zero Sum..." << endl;
    
    try {
        size_t l = 4;
        KZG::PublicKey pk = setup(3 * l);
        
        Fr omega = findPrimitiveRoot(l);
        
        Fr zero_sum = 0;
        vector<Fr> q = createValidSumCheckPolynomial(l, omega, zero_sum);
        
        bool result = sumCheck(pk, q, omega, l, zero_sum);
        
        if (result) {
            cout << "✓ Test 5 PASSED" << endl;
            return true;
        } else {
            cout << "✗ Test 5 FAILED" << endl;
            return false;
        }
    } catch (const exception& e) {
        cout << "✗ Test 5 FAILED with exception: " << e.what() << endl;
        return false;
    }
}

// Test case 6: Purely vanishing polynomial (sum should be 0)
bool testPureVanishingPolynomial() {
    cout << "Test 6: Pure Vanishing Polynomial..." << endl;
    
    try {
        size_t l = 4;
        KZG::PublicKey pk = setup(3 * l);
        
        Fr omega = findPrimitiveRoot(l);
        
        // Create the vanishing polynomial zh(x) = x^l - 1
        vector<Fr> q(l + 1, 0);
        q[0] = -1;  // -1
        q[l] = 1;   // x^l
        
        Fr expected_sum = 0;  // Should sum to 0 since it vanishes on H
        
        bool result = sumCheck(pk, q, omega, l, expected_sum);
        
        if (result) {
            cout << "✓ Test 6 PASSED" << endl;
            return true;
        } else {
            cout << "✗ Test 6 FAILED" << endl;
            return false;
        }
    } catch (const exception& e) {
        cout << "✗ Test 6 FAILED with exception: " << e.what() << endl;
        return false;
    }
}

// Test case 7: Verify NTT root properties
bool testNTTRootProperties() {
    cout << "Test 7: NTT Root Properties..." << endl;
    
    try {
        size_t l = 8;
        Fr omega = findPrimitiveRoot(l);
        
        // Verify omega^l = 1
        Fr omega_l;
        Fr::pow(omega_l, omega, l);
        Fr one = 1;
        
        if (omega_l != one) {
            cout << "✗ Test 7 FAILED: omega^l != 1" << endl;
            return false;
        }
        
        // Verify omega^(l/2) = -1 (for l > 2)
        if (l > 2) {
            Fr omega_half;
            Fr::pow(omega_half, omega, l / 2);
            Fr minus_one = -1;
            
            if (omega_half != minus_one) {
                cout << "✗ Test 7 FAILED: omega^(l/2) != -1" << endl;
                return false;
            }
        }
        
        // Verify all powers are distinct
        vector<Fr> powers(l);
        Fr curr = 1;
        for (size_t i = 0; i < l; i++) {
            powers[i] = curr;
            Fr::mul(curr, curr, omega);
        }
        
        for (size_t i = 0; i < l; i++) {
            for (size_t j = i + 1; j < l; j++) {
                if (powers[i] == powers[j]) {
                    cout << "✗ Test 7 FAILED: powers not distinct" << endl;
                    return false;
                }
            }
        }
        
        cout << "✓ Test 7 PASSED (NTT root properties verified)" << endl;
        return true;
        
    } catch (const exception& e) {
        cout << "✗ Test 7 FAILED with exception: " << e.what() << endl;
        return false;
    }
}

// Test case 8: Large domain stress test
bool testLargeDomainStress() {
    cout << "Test 8: Large Domain Stress Test (l=16)..." << endl;
    
    try {
        size_t l = 16;
        KZG::PublicKey pk = setup(4 * l);  // Even larger setup
        
        Fr omega = findPrimitiveRoot(l);
        
        Fr target_sum = 12345;
        vector<Fr> q = createValidSumCheckPolynomial(l, omega, target_sum);
        
        bool result = sumCheck(pk, q, omega, l, target_sum);
        
        if (result) {
            cout << "✓ Test 8 PASSED" << endl;
            return true;
        } else {
            cout << "✗ Test 8 FAILED" << endl;
            return false;
        }
    } catch (const exception& e) {
        cout << "✗ Test 8 FAILED with exception: " << e.what() << endl;
        return false;
    }
}

int main() {
    cout << "=== Univariate Sum Check PIOP Tests (with NTT) ===" << endl;
    cout << "Using BN_SNARK1 curve with MCL library and NTT" << endl << endl;
    
    // Initialize the curve
    initPairing(BN_SNARK1);
    
    int passed = 0;
    int total = 8;
    
    // Run all tests
    if (testBasicValidSumCheck()) passed++;
    cout << endl;
    
    if (testInvalidSum()) passed++;
    cout << endl;
    
    if (testNonVanishingPolynomial()) passed++;
    cout << endl;
    
    if (testLargerDomain()) passed++;
    cout << endl;
    
    if (testZeroSum()) passed++;
    cout << endl;
    
    if (testPureVanishingPolynomial()) passed++;
    cout << endl;
    
    if (testNTTRootProperties()) passed++;
    cout << endl;
    
    if (testLargeDomainStress()) passed++;
    cout << endl;
    
    // Summary
    cout << "=== Test Summary ===" << endl;
    cout << "Passed: " << passed << "/" << total << " tests" << endl;
    
    if (passed == total) {
        cout << "🎉 All tests PASSED!" << endl;
        return 0;
    } else {
        cout << "❌ Some tests FAILED!" << endl;
        return 1;
    }
}
