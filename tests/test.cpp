#include "zerotest.h"
#include "kzg.h"
#include "ntt.h"
#include <mcl/bn.hpp>
#include <iostream>
#include <vector>
#include <cassert>
#include <stdexcept>

using namespace std;
using namespace mcl;
using namespace bn;

// Test helper functions
void printTestResult(const string& testName, bool passed) {
    cout << "[" << (passed ? "PASS" : "FAIL") << "] " << testName << endl;
}

// Helper function to create a polynomial that vanishes on H = {1, w, w^2, ..., w^(l-1)}
vector<Fr> createVanishingPolynomial(Fr w, size_t l) {
    // Create polynomial q(x) = (x - 1)(x - w)(x - w^2)...(x - w^(l-1))
    vector<Fr> q = {-1, 1}; // Start with (x - 1)
    
    Fr curr = w;
    for (size_t i = 1; i < l; i++) {
        // Multiply by (x - w^i)
        vector<Fr> factor = {-curr, 1};
        vector<Fr> newQ(q.size() + 1, 0);
        
        for (size_t j = 0; j < q.size(); j++) {
            for (size_t k = 0; k < factor.size(); k++) {
                newQ[j + k] = newQ[j + k] + q[j] * factor[k];
            }
        }
        q = newQ;
        curr *= w;
    }
    
    return q;
}

// Helper function to create a random polynomial of given degree
vector<Fr> createRandomPolynomial(size_t degree) {
    vector<Fr> poly(degree + 1);
    for (size_t i = 0; i <= degree; i++) {
        poly[i].setByCSPRNG();
    }
    return poly;
}

// Test polynomial division function
bool testPolynomialDivision() {
    try {
        Fr w = findPrimitiveRoot(8);
        
        // Create a simple test case: divide x^2 - 1 by x - 1
        vector<Fr> dividend = {-1, 0, 1}; // x^2 - 1
        vector<Fr> divisor = {-1, 1};     // x - 1
        
        vector<Fr> quotient = polynomialDivision(dividend, divisor);
        
        // Expected result should be x + 1
        // Verify by checking that quotient * divisor = dividend
        bool correct = true;
        
        // For this simple case, we can manually verify
        // (x + 1)(x - 1) = x^2 - 1
        if (quotient.size() >= 2) {
            Fr expected_const = 1;  // constant term of x + 1
            Fr expected_linear = 1; // linear term of x + 1
            
            // Check if quotient represents x + 1 (allowing for some padding zeros)
            correct = (quotient[0] == expected_const && quotient[1] == expected_linear);
        }
        
        return correct;
    } catch (const exception& e) {
        cout << "Exception in testPolynomialDivision: " << e.what() << endl;
        return false;
    }
}

// Test zero test with a polynomial that should pass
bool testZeroTestValidCase() {
    try {
        // Setup KZG
        KZG::PublicKey pk = setup(16); // Setup with degree 16
        
        size_t l = 4;
        Fr w = findPrimitiveRoot(l);
        
        // Create a polynomial that vanishes on H
        vector<Fr> q = createVanishingPolynomial(w, l);
        
        // Test should return true
        bool result = zeroTest(pk, q, w, l);
        return result;
        
    } catch (const exception& e) {
        cout << "Exception in testZeroTestValidCase: " << e.what() << endl;
        return false;
    }
}

// Test zero test with a polynomial that should fail
bool testZeroTestInvalidCase() {
    try {
        // Setup KZG
        KZG::PublicKey pk = setup(16); // Setup with degree 16
        
        size_t l = 4;
        Fr w = findPrimitiveRoot(l);
        
        // Create a polynomial that does NOT vanish on H
        vector<Fr> q = {1, 2, 3}; // 3x^2 + 2x + 1
        
        // This should throw an exception since q doesn't vanish on H
        try {
            zeroTest(pk, q, w, l);
            return false; // Should not reach here
        } catch (const runtime_error& e) {
            // Expected exception
            return true;
        }
        
    } catch (const exception& e) {
        cout << "Unexpected exception in testZeroTestInvalidCase: " << e.what() << endl;
        return false;
    }
}

// Test with larger polynomial and domain
bool testZeroTestLargerDomain() {
    try {
        // Setup KZG
        KZG::PublicKey pk = setup(32); // Setup with larger degree
        
        size_t l = 8;
        Fr w = findPrimitiveRoot(l);
        
        // Create a polynomial that vanishes on H
        vector<Fr> q = createVanishingPolynomial(w, l);
        
        // Test should return true
        bool result = zeroTest(pk, q, w, l);
        return result;
        
    } catch (const exception& e) {
        cout << "Exception in testZeroTestLargerDomain: " << e.what() << endl;
        return false;
    }
}

// Test edge case with l = 1
bool testZeroTestSingleElement() {
    try {
        // Setup KZG
        KZG::PublicKey pk = setup(8); // Setup with sufficient degree
        
        size_t l = 1;
        Fr w = findPrimitiveRoot(2); // Need at least 2nd root of unity
        
        // For l=1, H = {1}, so q(x) should satisfy q(1) = 0
        vector<Fr> q = {-1, 1}; // x - 1
        
        bool result = zeroTest(pk, q, w, l);
        return result;
        
    } catch (const exception& e) {
        cout << "Exception in testZeroTestSingleElement: " << e.what() << endl;
        return false;
    }
}

// Test with polynomial that has extra factors beyond vanishing polynomial
bool testZeroTestWithExtraFactors() {
    try {
        // Setup KZG
        KZG::PublicKey pk = setup(32); // Setup with larger degree for complex polynomial
        
        size_t l = 4;
        Fr w = findPrimitiveRoot(l);
        
        // Create vanishing polynomial
        vector<Fr> vanishing = createVanishingPolynomial(w, l);
        
        // Multiply by an extra factor (x + 2)
        vector<Fr> extraFactor = {2, 1}; // x + 2
        vector<Fr> q(vanishing.size() + extraFactor.size() - 1, 0);
        
        for (size_t i = 0; i < vanishing.size(); i++) {
            for (size_t j = 0; j < extraFactor.size(); j++) {
                q[i + j] = q[i + j] + vanishing[i] * extraFactor[j];
            }
        }
        
        // This should still pass since q vanishes on H
        bool result = zeroTest(pk, q, w, l);
        return result;
        
    } catch (const exception& e) {
        cout << "Exception in testZeroTestWithExtraFactors: " << e.what() << endl;
        return false;
    }
}

int main() {
    cout << "Running Zero Test Implementation Tests..." << endl;
    cout << "=========================================" << endl;
    
    // Initialize MCL library
    initPairing(mcl::BN_SNARK1);
    
    int passedTests = 0;
    int totalTests = 0;
    
    // Run polynomial division tests
    totalTests++;
    bool divisionTest = testPolynomialDivision();
    printTestResult("Polynomial Division", divisionTest);
    if (divisionTest) passedTests++;
    
    // Run zero test with valid case
    totalTests++;
    bool validTest = testZeroTestValidCase();
    printTestResult("Zero Test - Valid Case", validTest);
    if (validTest) passedTests++;
    
    // Run zero test with invalid case
    // totalTests++;
    // bool invalidTest = testZeroTestInvalidCase();
    // printTestResult("Zero Test - Invalid Case (Should Throw)", invalidTest);
    // if (invalidTest) passedTests++;
    
    // Run zero test with larger domain
    totalTests++;
    bool largerTest = testZeroTestLargerDomain();
    printTestResult("Zero Test - Larger Domain", largerTest);
    if (largerTest) passedTests++;
    
    // Run zero test with single element
    totalTests++;
    bool singleTest = testZeroTestSingleElement();
    printTestResult("Zero Test - Single Element Domain", singleTest);
    if (singleTest) passedTests++;
    
    // Run zero test with extra factors
    totalTests++;
    bool extraFactorTest = testZeroTestWithExtraFactors();
    printTestResult("Zero Test - With Extra Factors", extraFactorTest);
    if (extraFactorTest) passedTests++;
    
    cout << "=========================================" << endl;
    cout << "Tests passed: " << passedTests << "/" << totalTests << endl;
    
    if (passedTests == totalTests) {
        cout << "All tests passed!" << endl;
        return 0;
    } else {
        cout << "Some tests failed!" << endl;
        return 1;
    }
}
