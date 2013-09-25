/*
 * =====================================================================================
 *
 *       Filename:  TestrationalNumber.cxx
 *
 *    Description:  Test for Rational Number class
 *
 *        Version:  1.0
 *        Created:  12/06/2012 04:51:16 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  Univesity of California, Berkeley
 *
 * =====================================================================================
 */

/* ------------------------------------------ 
 * Note on using operators for RationalNumber
 * ------------------------------------------
 *
 * Unfortunately overloaded operaotrs does not have appropriate proirity.
 * So use parenthesis as much.
 *
 * Wrong:   a *  b^c
 * Correct: c * (b^c)
 *
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
    RationalNumber a0;
    RationalNumber a1(2);
    RationalNumber a2(23.4564);
    RationalNumber a3(2,3);
    RationalNumber a4(16,-10);

    // Test Constructors
    std::cout << "Test Constructor" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout << "a0: " << a0 << std::endl;
    std::cout << "a1: (" << a1.GetNumerator() << "/" << a1.GetDenominator() << ")" << std::endl;
    std::cout << "a2: " << a2.GetAsString() << std::endl;
    std::cout << "a3: ";
    a3.Print();
    std::cout << "a4: ";
    RationalNumber::Print(std::cout,a4);
    std::cout << std::endl;

    // Test Accessors,Mutators
    std::cout << "Test Accessors, Mutators" << std::endl;
    std::cout << "------------------------" << std::endl;
    a3.SetNumerator(21);
    a3.SetDenominator(-9);
    std::cout << "a3: (" << a3.GetNumerator() << "/" << a3.GetDenominator() << ")" << std::endl;
    RationalNumber a31 = a3.NormalForm();
    std::cout << "a3 in normal form: " << a3.GetAsString() << std::endl;
    std::cout << std::endl;

    // Test Member Methods
    std::cout << "Test member methods" << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "a3 integer part: " << a3.IntegerPart() << std::endl;
    RationalNumber a32 = a3.FractionalPart();
    std::cout << "a3 fractional part: " << a32.GetAsString() << std::endl;
    RationalNumber a33 = a3.Inverse();
    std::cout << "a3 inverse: " << a33.GetAsString() << std::endl;
    RationalNumber a34 = a3.Opposite();
    std::cout << "a3 opposite: " << a34.GetAsString() << std::endl;
    std::cout << "GDC(6,-9): " << RationalNumber::GreatestCommonDivisor(-6,9) << std::endl;
    std::cout << "LCM(6,-9): " << RationalNumber::LeastCommonMultiple(6,-9) << std::endl;
    RationalNumber a5;
    a5.SetNumerator(21);
    a5.SetDenominator(-9);
    RationalNumber a51 = a5.Simplify();
    std::cout << "a5 simplify:: " << a51 << std::endl;
    RationalNumber a52 = a5.NormalForm();
    std::cout << "a5 normal form: " << a52 << std::endl;
    std::cout << "Evaluate a5: " << a5.Evaluate() << std::endl;
    std::cout << "Print a5: "; a5.Print();
    std::cout << "static Print a5: "; RationalNumber::Print(std::cout,a5);
    std::cout << std::endl;

    // Test Operators + - * / ^
    std::cout << "Test Operators" << std::endl;
    std::cout << "--------------" << std::endl;
    RationalNumber b1(7,6);
    RationalNumber b2(-4,9);
    RationalNumber b3 = b1 + b2;
    RationalNumber b4 = b1 - b2;
    RationalNumber b5 = b1 * b2;
    RationalNumber b6 = b1 / b2;
    RationalNumber b7 = b1^3;
    std::cout << "b1   : " << b1.GetAsString() << std::endl;
    std::cout << "b2   : " << b2.GetAsString() << std::endl;
    std::cout << "b1+b2: " << b3.GetAsString() << std::endl;
    std::cout << "b1-b2: " << b4.GetAsString() << std::endl;
    std::cout << "b1*b2: " << b5.GetAsString() << std::endl;
    std::cout << "b1/b2: " << b6.GetAsString() << std::endl;
    std::cout << "b1^3 : " << b7.GetAsString() << std::endl;
    std::cout << std::endl;

    // test operator =
    RationalNumber c1,c2;
    c2 = c1 = b1;
    std::cout << "c2 = c1 = b1: " << c1.GetAsString() << " and " << c2.GetAsString() << std::endl;
    
    // test operators += -= *= /=
    RationalNumber c3 = b1 += b2; std::cout << "b1 += b2: " << c3 << ", " << b1 << std::endl;
    RationalNumber c4 = b1 -= b2; std::cout << "b1 -= b2: " << c4 << ", " << b1 << std::endl;
    RationalNumber c5 = b1 *= b2; std::cout << "b1 *= b2: " << c5 << ", " << b1 << std::endl;
    RationalNumber c6 = b1 /= b2; std::cout << "b1 /= b2: " << c6 << ", " << b1 << std::endl;
    std::cout << std::endl;

    // test operaotrs < > <= >= == !=
    bool result;
    result = b1 >  b2;  std::cout << "b1 >  b2: " << result << std::endl;
    result = b1 <  b2;  std::cout << "b1 <  b2: " << result << std::endl;
    result = b1 <= b2;  std::cout << "b1 <= b2: " << result << std::endl;
    result = b1 <= b1;  std::cout << "b1 <= b1: " << result << std::endl;
    result = b2 >= b1;  std::cout << "b2 >= b1: " << result << std::endl;
    result = b2 >= b2;  std::cout << "b2 >= b2: " << result << std::endl;
    result = b1 == b2;  std::cout << "b1 == b2: " << result << std::endl;
    result = b1 == b1;  std::cout << "b1 == b1: " << result << std::endl;
    result = b1 != b1;  std::cout << "b1 != b1: " << result << std::endl;
    result = b1 != b2;  std::cout << "b1 != b2: " << result << std::endl;
    std::cout << std::endl;

    // test conversion operator
    int int_value = b1;
    double double_value = b1;
    std::cout << "int_value    = b1: " << int_value    << std::endl;
    std::cout << "double_value = b1: " << double_value << std::endl;
    std::cout << std::endl;

    // Test templates for operators
    std::cout << "Test template for operators" << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::cout << "b1      : " << b1       << std::endl;
    std::cout << "b1 + 2  : " << b1 + 2   << std::endl;
    std::cout << "b1 + 2.1: " << b1 + 2.1 << std::endl;
    std::cout << "b1 - 2  : " << b1 - 2   << std::endl;
    std::cout << "b1 - 2.1: " << b1 - 2.1 << std::endl;
    std::cout << "b1 * 2  : " << b1 + 2.1 << std::endl;
    std::cout << "b1 * 2.1: " << b1 * 2   << std::endl;
    std::cout << "b1 / 2  : " << b1 / 2   << std::endl;
    std::cout << "b1 / 2.1: " << b1 / 2.1 << std::endl;
    std::cout << std::endl;

    c1 = c2 = 2;   std::cout << "c1 = c2 = 2: "   << c1 << ", " << c2 << std::endl;
    c1 = c2 = 2.1; std::cout << "c1 = c2 = 2.1: " << c1 << ", " << c2 << std::endl;
    c1 += 2;       std::cout << "c1 += 2  : "     << c1 << std::endl;
    c1 += 2.1;     std::cout << "c1 += 2.1: "     << c1 << std::endl;
    c1 -= 2.1;     std::cout << "c1 -= 2.1: "     << c1 << std::endl;
    c1 -= 2;       std::cout << "c1 -= 2  : "     << c1 << std::endl;
    c1 *= 2;       std::cout << "c1 *= 2  : "     << c1 << std::endl;
    c1 *= 2.1;     std::cout << "c1 *= 2.1: "     << c1 << std::endl;
    c1 /= 2;       std::cout << "c1 /= 2  : "     << c1 << std::endl;
    c1 /= 2.1;     std::cout << "c1 /= 2.1: "     << c1 << std::endl;
    std::cout << std::endl;

    c1 = -2.3;           std::cout << "c1: "          << c1     << std::endl;
    result = c1 <   2.1; std::cout << "c1 <   2.1: "  << result << std::endl;
    result = c1 >   2.1; std::cout << "c1 >   2.1: "  << result << std::endl;
    result = c1 <=  2.1; std::cout << "c1 <=  2.1: "  << result << std::endl;
    result = c1 <= -2.3; std::cout << "c1 <= -2.3: "  << result << std::endl;
    result = c1 >=  2.1; std::cout << "c1 >=  2.1: "  << result << std::endl;
    result = c1 >= -2.3; std::cout << "c1 >= -2.3: "  << result << std::endl;
    result = c1 == -2.3; std::cout << "c1 == -2.3: "  << result << std::endl;
    result = c1 !=  2.1; std::cout << "c1 !=  2.1: "  << result << std::endl;
    result = c1 != -2.3; std::cout << "c1 != -2.3: "  << result << std::endl;
    std::cout << std::endl;

    // Test Left Hand side Operators
    std::cout << "Test left hand side opertors" << std::endl;
    std::cout << "----------------------------" << std::endl;
    std::cout << "2   + c1: " << 2   + c1 << std::endl;
    std::cout << "2.1 + c1: " << 2.1 + c1 << std::endl;
    std::cout << "2   - c1: " << 2   - c1 << std::endl;
    std::cout << "2.1 - c1: " << 2.1 - c1 << std::endl;
    std::cout << "2   * c1: " << 2   * c1 << std::endl;
    std::cout << "2.1 * c1: " << 2.1 * c1 << std::endl;
    std::cout << "2   / c1: " << 2   / c1 << std::endl;
    std::cout << "2.1 / c1: " << 2.1 / c1 << std::endl;
    std::cout << std::endl;

    result = 2    < c1;  std::cout << "2    <  c1: " << result << std::endl;
    result = 2.1  < c1;  std::cout << "2.1  <  c1: " << result << std::endl;
    result = 2    > c1;  std::cout << "2    >  c1: " << result << std::endl;
    result = 2.1  > c1;  std::cout << "2.1  >  c1: " << result << std::endl;
    result = 2    <= c1; std::cout << "2    <= c1: " << result << std::endl;
    result = -2.3 <= c1; std::cout << "-2.3 <= c1: " << result << std::endl;
    result = 2    >= c1; std::cout << "2    >= c1: " << result << std::endl;
    result = -2.3 >= c1; std::cout << "-2.3 >= c1: " << result << std::endl;
    result = 2    == c1; std::cout << "2    == c1: " << result << std::endl;
    result = -2.3 == c1; std::cout << "-2.3 == c1: " << result << std::endl;
    result = 2    != c1; std::cout << "2    != c1: " << result << std::endl;
    result = -2.3 != c1; std::cout << "-2.3 != c1: " << result << std::endl; 
    std::cout << std::endl;

    return 0;
}
