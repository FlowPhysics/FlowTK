/*
 * =====================================================================================
 *
 *       Filename:  TestPolynomial.cxx
 *
 *    Description:  Test for Polynomial class
 *
 *        Version:  1.0
 *        Created:  12/20/2012 08:32:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

/* ---------------------------------------
 * Note on using operators for Polynomials
 * ---------------------------------------
 *
 * Polynomial operators with RHS RationalNumbers is being mistaken with
 * RationalNumber operators with LHS template polynomial
 *
 * asusme P is Polynomial and RN is RationalNumber. We want P + RN
 *
 * 1- P.operator+(RN)  this is templated RHS operator for Polynomial
 * 2- operator+(P,RN)  this is templated LHS operator friend of RationalNumber
 *
 * To avoid this always use the following:
 * P + Polynomial(RN)  now it is non-templated RHS operator for Polynomial + Polynomial
 */

// =======
// Headers
// =======

#include <GeneralMath.h>
#include <iostream>

// ====
// Main
// ====

int main(int argc, char *argv[])
{
    // Test Constructor
    std::cout << "Test Constructor" << std::endl;
    std::cout << "----------------" << std::endl;
    Polynomial P0;
    std::cout << "P0 Degree: " << P0.GetDegree() << std::endl;
    std::cout << "P0: ";
    Polynomial::Print(std::cout,P0);

    Polynomial P1(-2);
    std::cout << "P1 Degree: " << P1.GetDegree() << std::endl;
    std::cout << "P1: " << P1 << std::endl;

    Polynomial P2(1.2);
    std::cout << "P2 Degree: " << P2.GetDegree() << std::endl;
    for(unsigned int i=0; i<=P2.GetDegree(); i++)
    {
        std::cout << "P2 Coefficient[" << i << "]: " << P2.GetCoefficient(i).GetAsString() << std::endl;
    }

    RationalNumber r(3,4);
    Polynomial P3(r);
    std::cout << "P3 Degree: " << P3.GetDegree() << std::endl;
    std::cout << "P3: ";
    P3.Print();
    std::cout << std::endl;

    // Test SetCoefficient
    Polynomial P4;
    RationalNumber r1(3,4);
    RationalNumber r2(10,4);
    P4.SetDegree(4);
    P4.SetCoefficient(0,2);
    P4.SetCoefficient(1,4);
    P4.SetCoefficient(2,r1);
    P4.SetCoefficient(3,r2);
    P4.SetCoefficient(4,RationalNumber(-15,9));
    std::cout << "P4 Degree: " << P4.GetDegree() << std::endl;
    std::cout << "P4: ";
    P4.Print();

    // Test Copy Constructor
    Polynomial P5(P4);
    std::cout << "P5 Degree: " << P5.GetDegree() << std::endl;
    std::cout << "P5: ";
    Polynomial::Print(std::cout,P5);
    std::cout << std::endl;

    // Test Accessors
    std::cout << "Test Accessors" << std::endl;
    std::cout << "--------------" << std::endl;
    P5.SetCoefficient(1,1.6);
    P5.AddNextCoefficient(3);
    P5.AddNextCoefficient(-2.4);
    P5.AddNextCoefficient(RationalNumber(-6,4));
    P5.SetCoefficient(2,RationalNumber(7,3));
    std::cout << "P5; " << P5 << std::endl;
    P5.SetCoefficient(10,RationalNumber(4,5));
    std::cout << "P5: " << P5 << std::endl;
    std::cout << std::endl;
    
    // Test change Degree
    std::cout << "Test change Degree" << std::endl;
    std::cout << "------------------" << std::endl;
    P4.SetDegree(7);
    std::cout << "P4: " << P4 << std::endl;
    P4.SetDegree(5);
    std::cout << "P4: " << P4 << std::endl;
    std::cout << std::endl;

    // Test Differentiate
    std::cout << "Test Differentiate" << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "P4       : " << P4 << std::endl;
    std::cout << "Diff-0 P4: " << P4.Differentiate(0) << std::endl;
    Polynomial D1_P4 = P4.Differentiate();
    std::cout << "Diff-1 P4: " << D1_P4 << std::endl; 
    Polynomial D2_P4 = P4.Differentiate(2);
    std::cout << "Diff-2 P4: " << D2_P4 << std::endl;
    Polynomial D3_P4 = Polynomial::Differentiate(P4,3);
    std::cout << "Diff-3 P4: " << D3_P4 << std::endl;
    Polynomial D4_P4 = Polynomial::Differentiate(D3_P4);
    std::cout << "Diff-4 P4: " << D4_P4 << std::endl;
    Polynomial D5_P4 = Polynomial::Differentiate(P4,5);
    std::cout << "D5_P4    : " << D5_P4 << std::endl;
    std::cout << std::endl;

    // Test Integrate
    std::cout << "Test Integrate" << std::endl;
    std::cout << "--------------" << std::endl;
    std::cout << "P4    : " << P4 << std::endl;
    std::cout << "I-0 P4: " << P4.Integrate(0) << std::endl;
    Polynomial I1_P4 = P4.Integrate();
    std::cout << "I-1 P4: " << I1_P4 << std::endl;
    std::cout << "I-2 P4: " << Polynomial::Integrate(I1_P4) << std::endl;
    std::cout << "I-3 P4: " << P4.Integrate(3) << std::endl;
    std::cout << "I-4 P4: " << Polynomial::Integrate(P4,4) << std::endl;

    // Test Evaluate
    std::cout << "Test Evaluate" << std::endl;
    std::cout << "-------------" << std::endl;
    std::cout << "P4: " << P4 << std::endl;
    std::cout << P4.Evaluate(2) << std::endl;
    std::cout << Polynomial::Evaluate(P4,-2.4) << std::endl;
    std::cout << std::endl;

    // Test Resize
    std::cout << "Test Resize" << std::endl;
    std::cout << "-----------" << std::endl;
    std::cout << "P4        : " << P4 << std::endl;
    P4.Resize(7);
    std::cout << "resize  P4: " << P4 << std::endl;
    P4.Resize(6);
    std::cout << "resize  P4: " << P4 << std::endl;
    P4.Extend(8);
    std::cout << "extend  P4: " << P4 << std::endl;
    P4.Squeeze(5);
    std::cout << "squeeze P4: " << P4 << std::endl;
    std::cout << std::endl;

    // test RHS operators
    P5.Resize(3);
    std::cout << "P4     : " << P4 << std::endl;
    std::cout << "P5     : " << P5 << std::endl;
    std::cout << "P4 + P5: " << P4 + P5 << std::endl;
    std::cout << "P4 - P5: " << P4 - P5 << std::endl;
    std::cout << "P4 * P5: " << P4 * P5 << std::endl;
    Polynomial P6 = P5 ^ 3;
    std::cout << "P5 ^ 3 : " << P6  << std::endl;
    std::cout << std::endl;

    Polynomial P7;
    P6 = P7 = P5;
    std::cout << "P5           : " << P5 << std::endl;
    std::cout << "P7 = P5      : " << P7 << std::endl;
    std::cout << "P6 = P7 = P5 : " << P6 << std::endl;
    std::cout << std::endl;

    P6 += P4; std::cout << "P6 += P4: " << P6 << std::endl;
    P6 -= P4; std::cout << "P6 -= P4: " << P6 << std::endl;
    P6 *= P4; std::cout << "P6 *= P4: " << P6 << std::endl;
    std::cout << "P6[3]:  " << P6[3] <<  std::endl;
    std::cout << std::endl;

    bool result;
    result = P4 == P5; std::cout << "P4 == P5: " << result << std::endl;
    result = P4 == P4; std::cout << "P4 == P4: " << result << std::endl;
    result = P4 != P5; std::cout << "P4 != P5: " << result << std::endl;
    result = P4 != P4; std::cout << "P4 != P4: " << result << std::endl;
    std::cout << std::endl;

    P6.SetCoefficient(0,-3.6);
    std::cout << "P6: " << P6 << std::endl;
    int integer_value;
    double double_value;
    integer_value = P6; std::cout << "integer_value = P6: " << integer_value << std::endl;
    double_value  = P6; std::cout << "double_value  = P6: " << double_value  << std::endl;
    std::cout << std::endl;

    // Test Template RHS Operators
    std::cout << "Test template RHS operators" << std::endl; 
    std::cout << "---------------------------" << std::endl;
    std::cout << "P4       : "   << P4       << std::endl;
    std::cout << "P4 + 2   : "   << P4 + 2   << std::endl;
    std::cout << "P4 + 2.1 : "   << P4 + 2.1 << std::endl;
    std::cout << "P4 - 2   : "   << P4 - 2   << std::endl;
    std::cout << "P4 - 2.1 : "   << P4 - 2.1 << std::endl;
    std::cout << "P4 * 2   : "   << P4 * 2   << std::endl;
    std::cout << "P4 * 2.1 : "   << P4 * 2.1 << std::endl;
    std::cout << std::endl;

    P1 = 2;                    std::cout << "P1 = 2     : " << P1 << std::endl;
    P1 = 2.1;                  std::cout << "P1 = 2.1   ; " << P1 << std::endl;
    P1 = RationalNumber(2,3);  std::cout << "P1 = (2/3) : " << P1 << std::endl;
    P4 += 2;                   std::cout << "P4 += 2    : " << P4 << std::endl;
    P4 += 2.1;                 std::cout << "P4 += 2.1  : " << P4 << std::endl;
    P4 += RationalNumber(2,3); std::cout << "P4 += (2/3): " << P4 << std::endl;
    P4 -= 2;                   std::cout << "P4 -= 2    : " << P4 << std::endl;
    P4 -= 2.1;                 std::cout << "P4 -= 2.1  : " << P4 << std::endl;
    P4 -= RationalNumber(2,3); std::cout << "P4 -= (2/3): " << P4 << std::endl;
    P4 *= 2;                   std::cout << "P4 *= 2    : " << P4 << std::endl;
    P4 *= 2.1;                 std::cout << "P4 *= 2.1  : " << P4 << std::endl;
    P4 *= RationalNumber(2,3); std::cout << "P4 *= (2/3): " << P4 << std::endl;
    std::cout << std::endl;

    result = P4 == 2;                   std::cout << "P4 == 2    : " << result << std::endl;
    result = P4 == 2.1;                 std::cout << "P4 == 2.1  : " << result << std::endl;
    result = P4 == RationalNumber(2,3); std::cout << "P4 == (2/3): " << result << std::endl;
    result = P4 != 2;                   std::cout << "P4 != 2    : " << result << std::endl;
    result = P4 != 2.1;                 std::cout << "P4 != 2.1  : " << result << std::endl;
    result = P4 != RationalNumber(2,3); std::cout << "P4 != (2/3): " << result << std::endl;
    std::cout << std::endl;

    // Test Template LHS Operators
    std::cout << "P5        : " << P5 << std::endl;
    std::cout << "2     + P5: " << 2                   + P5 << std::endl; 
    std::cout << "2.1   + P5: " << 2.1                 + P5 << std::endl; 
    std::cout << "(2/3) + P5: " << RationalNumber(2,3) + P5 << std::endl; 
    std::cout << "2     - P5: " << 2                   - P5 << std::endl; 
    std::cout << "2.1   - P5: " << 2.1                 - P5 << std::endl; 
    std::cout << "(2/3) - P5: " << RationalNumber(2,3) - P5 << std::endl; 
    std::cout << "2     * P5: " << 2                   * P5 << std::endl; 
    std::cout << "2.1   * P5: " << 2.1                 * P5 << std::endl; 
    std::cout << "(2/3) * P5: " << RationalNumber(2,3) * P5 << std::endl; 
    std::cout << std::endl;

    result = 2                   == P5; std::cout << "2     == P5 : " << result << std::endl;
    result = 2.1                 == P5; std::cout << "2.1   == P5 : " << result << std::endl;
    result = RationalNumber(2,3) == P5; std::cout << "(2/3) == P5 : " << result << std::endl;
    result = 2                   != P5; std::cout << "2     != P5 : " << result << std::endl;
    result = 2.1                 != P5; std::cout << "2.1   != P5 : " << result << std::endl;
    result = RationalNumber(2,3) != P5; std::cout << "(2/3) != P5 : " << result << std::endl;
    std::cout << std::endl;

    // PROBLEMS TO BE SOLVED
    // Polynomial RHS Operators with RationalNumbers is being mistaken with 
    // RationalNumber LHS operators with template.
    std::cout << "Operators with RationalNumber in RHS" << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "THIS IS WRONG:" << std::endl;
    std::cout << "--------------" << std::endl;
    std::cout << "P4 + (2/3): " << P4 + RationalNumber(2,3) << std::endl;
    std::cout << "P4 - (2/3): " << P4 - RationalNumber(2,3) << std::endl;
    std::cout << "P4 * (2/3): " << P4 * RationalNumber(2,3) << std::endl;
    std::cout << std::endl;

    std::cout << "THIS IS CORRECT:" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout << "P4 + (2/3): " << P4 + Polynomial(RationalNumber(2,3)) << std::endl;
    std::cout << "P4 - (2/3): " << P4 - Polynomial(RationalNumber(2,3)) << std::endl;
    std::cout << "P4 * (2/3): " << P4 * Polynomial(RationalNumber(2,3)) << std::endl;
    std::cout << std::endl;

    return 0;
}
