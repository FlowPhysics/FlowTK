/*
 * =====================================================================================
 *
 *       Filename:  GeneralMath.h
 *
 *    Description:  GeneralMath Library
 *
 *        Version:  1.0
 *        Created:  12/04/2012 02:55:20 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __GeneralMath_h
#define __GeneralMath_h

// ======
// Macros
// ======

#ifndef HERE
#define HERE std::cout << "DEBUG: " << __FILE__ << " at " << __LINE__ << std::endl;
#endif

// ====================
// Forward Declarations
// ====================

// Complete types
#include <iostream>

// Incomplete types
class Combinatorics;
class Vector;
class Matrix;
class Polynomial;
class RationalNumber;

// ============
// General Math
// ============

class GeneralMath
{
    public:
        GeneralMath();
        ~GeneralMath();
};

// =============
// Combinatorics
// =============

class Combinatorics
{
    public:
        Combinatorics();
        ~Combinatorics();

        // Methods
        static unsigned long int Factorial(unsigned int integer);
        static unsigned long int Permutation(unsigned int n, unsigned int k);
        static unsigned long int Combination(unsigned int n, unsigned int k);
};

// ======
// Vector
// ======

class Vector
{
    public:
        // Constructors
        Vector();
        ~Vector();

        // Methods
        static double FindMaximumValue(const double *VectorArray, unsigned int VectorArrayLength);

    private:
        RationalNumber *Components;
        unsigned int Length;
};

// ======
// Matrix
// ======

class Matrix
{
    public:
        Matrix();
        ~Matrix();

        // Methods
        static double DeterminantOf2x2Matrix(const double Matrix[2][2]);
        static double DeterminantOf2x2Matrix(const double Column1[2],const double Column2[2]);
        static double DeterminantOf2x2Matrix(
                const double & Matrix_11, const double & Matrix_12,
                const double & Matrix_21, const double & Matrix_22);

        static double DeterminantOf3x3Matrix(const double Matrix[3][3]);
        static double DeterminantOf3x3Matrix(const double Column1[3], const double Column2[3], const double Column3[3]);
        static double DeterminantOf3x3Matrix(
                const double & Matrix_11, const double & Matrix_12, const double & Matrix_13,
                const double & Matrix_21, const double & Matrix_22, const double & Matrix_23,
                const double & Matrix_31, const double & Matrix_32, const double & Matrix_33);

        static double DeterminantOfNxNMatrix(double **Matrix,unsigned int MatrixSize);
        static double DeterminantOfNxNMatrix(const double *VectorizedMatrix,unsigned int MatrixSize);
 
        static void DirectSolve2x2LinearEquations(
                const double CoefficientMatrix[2][2],
                const double KnownVector[2],
                double UnknownVector[2]);

        static void DirectSolve2x2TransposeLinearEquations(
                double TransposedMatrix[2][2],
                double KnownVector[2],
                double UnknownVector[2]);
        
        static void DirectSolve3x3LinearEquations(
                const double CoefficientMatrix[3][3],
                const double KnwonVector[3],
                double UnknownVector[3]);

        static void DirectSolve3x3TransposeLinearEquations(
                double TransposedMatrix[3][3],
                double KnownVector[3],
                double UnknownVector[3]);

        static void DirectSolveNxNLinearEquations(
                double **CoefficientMatrix,
                double *KnownVector,
                unsigned int NumberOfEquations,
                double *UnknownVectors);

        static void DirectSolveNxNTransposeLinearEquations(
                double **TransposedMatrix,
                double *KnownVector,
                unsigned int NumberOfEquations,
                double *UnknownVector);

        static void IterativeSolveNxNLinearEquations(
                const double **CoefficientMatrix,
                const double *KnownVector,
                unsigned int NumberOfEquations,
                double *UnknownVector);
};

// ==========
// Polynomial
// ==========

class Polynomial
{
    public:
        // Constructors
        Polynomial();
        Polynomial(int integer);
        Polynomial(double Double);
        Polynomial(RationalNumber rational);
        Polynomial(const Polynomial & rhs);
        ~Polynomial();

        // Accessors, Mutators
        RationalNumber * GetCoefficients() const;
        RationalNumber GetCoefficient(unsigned int index) const;
        void SetCoefficient(unsigned int index, double coef);
        void SetCoefficient(unsigned int index, RationalNumber RNcoef);
        void SetCoefficients(double * DoubleCoefs, unsigned int degree);
        void SetCoefficients(RationalNumber * RNCoefs, unsigned int degree);
        void AddNextCoefficient(double coef);
        void AddNextCoefficient(RationalNumber RNcoef);
        unsigned int GetDegree() const;
        void SetDegree(unsigned int degree);
        const char * GetAsString() const;

        // Methods
        Polynomial Differentiate();
        static Polynomial Differentiate(const Polynomial & P);
        Polynomial Differentiate(unsigned int order);
        static Polynomial Differentiate(const Polynomial & P,unsigned int order);
        Polynomial Integrate();
        static Polynomial Integrate(const Polynomial & P);
        Polynomial Integrate(unsigned int order);
        static Polynomial Integrate(const Polynomial & P, unsigned int order);
        RationalNumber Evaluate(const RationalNumber & rational);
        static RationalNumber Evaluate(const Polynomial & P, const RationalNumber & rational);
        void Resize(unsigned int NewDegree);
        void Squeeze();
        void Squeeze(unsigned int NewDegree);
        void Extend(unsigned int NewDegree);

        void Print() const;
        static void Print(std::ostream & os, const Polynomial & P);

        // RHS Operators
        Polynomial operator+(const Polynomial & rhs) const;
        Polynomial operator-(const Polynomial & rhs) const;
        Polynomial operator*(const Polynomial & rhs) const;
        Polynomial operator^(const int power) const;
        Polynomial & operator=(const Polynomial & rhs);
        Polynomial & operator+=(const Polynomial & rhs);
        Polynomial & operator-=(const Polynomial & rhs);
        Polynomial & operator*=(const Polynomial & rhds);
        RationalNumber operator[](unsigned int offset) const;
        bool operator==(const Polynomial & rhs) const;
        bool operator!=(const Polynomial & rhs) const;
        operator int() const;
        operator double() const;
        operator RationalNumber() const;

        // RHS Template operators
        template <class Type> Polynomial operator+(const Type & rhs) const;
        template <class Type> Polynomial operator-(const Type & rhs) const;
        template <class Type> Polynomial operator*(const Type & rhs) const;
        template <class Type> Polynomial & operator=(const Type & rhs);
        template <class Type> Polynomial & operator+=(const Type & rhs);
        template <class Type> Polynomial & operator-=(const Type & rhs);
        template <class Type> Polynomial & operator*=(const Type & rhs);
        template <class Type> bool operator==(const Type & rhs) const;
        template <class Type> bool operator!=(const Type & rhs) const;

        // LHS operators
        friend std::ostream & operator<<(std::ostream & os, const Polynomial & rhs);
        template <class Type> friend Polynomial operator+(const Type & lhs, const Polynomial & rhs);
        template <class Type> friend Polynomial operator-(const Type & lhs, const Polynomial & rhs);
        template <class Type> friend Polynomial operator*(const Type & lhs, const Polynomial & rhs);
        template <class Type> friend bool operator==(const Type & lhs, const Polynomial & rhs);
        template <class Type> friend bool operator!=(const Type & lhs, const Polynomial & rhs);

    protected:
        RationalNumber *Coefficients;
        unsigned short int Degree;
};

// =================================
// Polynomial RHS Operator Templates
// =================================

// +
template <class Type> Polynomial Polynomial::operator+(const Type & rhs) const
{
    return this->operator+(Polynomial(rhs));
}

// -
template <class Type> Polynomial Polynomial::operator-(const Type & rhs) const
{
    return this->operator-(Polynomial(rhs));
}

// *
template <class Type> Polynomial Polynomial::operator*(const Type & rhs) const
{
    return this->operator*(Polynomial(rhs));
}

// =
template <class Type> Polynomial & Polynomial::operator=(const Type & rhs)
{
    return this->operator=(Polynomial(rhs));
}

// +=
template <class Type> Polynomial & Polynomial::operator+=(const Type & rhs)
{
    return this->operator+=(Polynomial(rhs));
}

// -=
template <class Type> Polynomial & Polynomial::operator-=(const Type & rhs)
{
    return this->operator-=(Polynomial(rhs));
}

// *=
template <class Type> Polynomial & Polynomial::operator*=(const Type & rhs)
{
    return this->operator*=(Polynomial(rhs));
}

// ==
template <class Type> bool Polynomial::operator==(const Type & rhs) const
{
    return this->operator==(Polynomial(rhs));
}

// !=
template <class Type> bool Polynomial::operator!=(const Type & rhs) const
{
    return this->operator!=(Polynomial(rhs));
}

// =================================
// Polynomial LHS Template operators
// =================================

// +
template <class Type> Polynomial operator+(const Type & lhs, const Polynomial & rhs)
{
    return Polynomial(rhs) + (Polynomial(lhs));
}

// -
template <class Type> Polynomial operator-(const Type & lhs, const Polynomial & rhs)
{
    return (-1 * rhs) + Polynomial(lhs);
}

// *
template <class Type> Polynomial operator*(const Type & lhs, const Polynomial & rhs)
{
    return rhs * Polynomial(lhs);
}

// ==
template <class Type> bool operator==(const Type & lhs, const Polynomial & rhs)
{
    return rhs == Polynomial(lhs);
}

// !=
template <class Type> bool operator!=(const Type & lhs, const Polynomial & rhs)
{
    return rhs != Polynomial(lhs);
}

// ===============
// Rational Number
// ===============

class RationalNumber
{
    public:
        // Constructors
        RationalNumber();
        RationalNumber(double Double);
        RationalNumber(long int numerator, long int denominator);
        RationalNumber(const RationalNumber & rhs);
        ~RationalNumber();

        // Accessors, Mutators
        long int GetNumerator() const { return this->Numerator; };
        void SetNumerator(long int numerator) { this->Numerator = numerator; }

        long int GetDenominator() const { return this->Denominator; }
        void SetDenominator(long int denominator);

        const char * GetAsString() const;

        // Member Methods
        long int IntegerPart() const;
        RationalNumber FractionalPart() const;
        RationalNumber Inverse() const;
        RationalNumber Opposite() const;
        static long int GreatestCommonDivisor(long int a, long int b);
        static long int LeastCommonMultiple(long int a, long int b);
        RationalNumber & Simplify();
        RationalNumber &  NormalForm();

        double Evaluate() const;
        void Print() const;
        static void Print(std::ostream & os,const RationalNumber & RN);

        // RHS Operators
        RationalNumber operator+(const RationalNumber & rhs) const;
        RationalNumber operator-(const RationalNumber & rhs) const;
        RationalNumber operator*(const RationalNumber & rhs) const;
        RationalNumber operator/(const RationalNumber & rhs) const;
        RationalNumber operator^(int power) const;
        RationalNumber & operator=(const RationalNumber & rhs);
        RationalNumber & operator+=(const RationalNumber & rhs);
        RationalNumber & operator-=(const RationalNumber & rhs);
        RationalNumber & operator*=(const RationalNumber & rhs);
        RationalNumber & operator/=(const RationalNumber & rhs);
        bool operator<(const RationalNumber & rhs) const;
        bool operator>(const RationalNumber & rhs) const;
        bool operator<=(const RationalNumber & rhs) const;
        bool operator>=(const RationalNumber & rhs) const;
        bool operator==(const RationalNumber & rhs) const;
        bool operator!=(const RationalNumber & rhs) const;
        operator int() const;
        operator double() const;

        // RHS Template operators
        template <class Type> RationalNumber operator+(const Type & rhs) const;
        template <class Type> RationalNumber operator-(const Type & rhs) const;
        template <class Type> RationalNumber operator*(const Type & rhs) const;
        template <class Type> RationalNumber operator/(const Type & rhs) const;
        template <class Type> RationalNumber & operator=(const Type & rhs);
        template <class Type> RationalNumber & operator+=(const Type & rhs);
        template <class Type> RationalNumber & operator-=(const Type & rhs);
        template <class Type> RationalNumber & operator*=(const Type & rhs);
        template <class Type> RationalNumber & operator/=(const Type & rhs);
        template <class Type> bool operator<(const Type & rhs) const;
        template <class Type> bool operator>(const Type & rhs) const;
        template <class Type> bool operator<=(const Type & rhs) const;
        template <class Type> bool operator>=(const Type & rhs) const;
        template <class Type> bool operator==(const Type & rhs) const;
        template <class Type> bool operator!=(const Type & rhs) const;

        // LHS operators
        friend std::ostream & operator<<(std::ostream & os,const RationalNumber & rhs) ;
        template <class Type> friend RationalNumber operator+(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend RationalNumber operator-(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend RationalNumber operator*(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend RationalNumber operator/(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend bool operator<(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend bool operator>(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend bool operator<=(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend bool operator>=(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend bool operator==(const Type & lhs, const RationalNumber & rhs);
        template <class Type> friend bool operator!=(const Type & lhs, const RationalNumber & rhs);

    protected:
        // Member Methods
        bool CheckValid() const;
        void AdjustSign();

        // Member Data
        long int Numerator;
        long int Denominator;
};

// =====================================
// RationalNumber RHS Operator Templates
// =====================================

// +
template <class Type> RationalNumber RationalNumber::operator+(const Type & rhs) const
{
    return this->operator+(RationalNumber(rhs));
}

// -
template <class Type> RationalNumber RationalNumber::operator-(const Type & rhs) const
{
    return this->operator-(RationalNumber(rhs));
}

// *
template <class Type> RationalNumber RationalNumber::operator*(const Type & rhs) const
{
    return this->operator*(RationalNumber(rhs));
}

// /
template <class Type> RationalNumber RationalNumber::operator/(const Type & rhs) const
{
    return this->operator/(RationalNumber(rhs));
}

// =
template <class Type> RationalNumber & RationalNumber::operator=(const Type & rhs)
{
    return this->operator=(RationalNumber(rhs));
}

// +=
template <class Type> RationalNumber & RationalNumber::operator+=(const Type & rhs)
{
    return this->operator+=(RationalNumber(rhs));
}

// -=
template <class Type> RationalNumber & RationalNumber::operator-=(const Type & rhs)
{
    return this->operator-=(RationalNumber(rhs));
}

// *=
template <class Type> RationalNumber & RationalNumber::operator*=(const Type & rhs)
{
    return this->operator*=(RationalNumber(rhs));
}

// /=
template <class Type> RationalNumber & RationalNumber::operator/=(const Type & rhs)
{
    return this->operator/=(RationalNumber(rhs));
}

// <
template <class Type> bool RationalNumber::operator<(const Type & rhs) const
{
    return this->operator<(RationalNumber(rhs));
}

// >
template <class Type> bool RationalNumber::operator>(const Type & rhs) const
{
    return this->operator>(RationalNumber(rhs));
}

// <=
template <class Type> bool RationalNumber::operator<=(const Type & rhs) const
{
    return this->operator<=(RationalNumber(rhs));
}

// >=
template <class Type> bool RationalNumber::operator>=(const Type & rhs) const
{
    return this->operator>=(RationalNumber(rhs));
}

// ==
template <class Type> bool RationalNumber::operator==(const Type & rhs) const
{
    return this->operator==(RationalNumber(rhs));
}

// !=
template <class Type> bool RationalNumber::operator!=(const Type & rhs) const
{
    return this->operator!=(RationalNumber(rhs));
}

// ====================================
// RatonalNumner LHS Template Operators
// ====================================

// +
template <class Type> RationalNumber operator+(const Type & lhs, const RationalNumber & rhs)
{
    return rhs + lhs;
}

// -
template <class Type> RationalNumber operator-(const Type & lhs, const RationalNumber & rhs)
{
    return (-1 * rhs) + lhs;
}

// *
template <class Type> RationalNumber operator*(const Type & lhs, const RationalNumber & rhs)
{
    return rhs * lhs;
}

// /
template <class Type> RationalNumber operator/(const Type & lhs, const RationalNumber & rhs)
{
    return lhs * rhs.Inverse();
}

// <
template <class Type> bool operator<(const Type & lhs, const RationalNumber & rhs)
{
    return rhs > lhs;
}

// >
template <class Type> bool operator>(const Type & lhs, const RationalNumber & rhs)
{
    return rhs < lhs;
}

// <=
template <class Type> bool operator<=(const Type & lhs, const RationalNumber & rhs)
{
    return rhs >= lhs;
}

// >=
template <class Type> bool operator>=(const Type & lhs, const RationalNumber & rhs)
{
    return rhs <= lhs;
}

// ==
template <class Type> bool operator==(const Type & lhs, const RationalNumber & rhs)
{
    return rhs == lhs;
}

// !=
template <class Type> bool operator!=(const Type & lhs, const RationalNumber & rhs)
{
    return rhs != lhs;
}

#endif
