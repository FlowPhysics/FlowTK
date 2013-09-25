/*
 * =====================================================================================
 *
 *       Filename:  GeneralMath.cxx
 *
 *    Description:  GeneralMath Library
 *
 *        Version:  1.0
 *        Created:  12/04/2012 02:54:42 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

// =======
// Headers
// =======

#include <GeneralMath.h>
#include <iostream>
#include <math.h>
#include <sstream>

// ======
// Macros
// ======

#define TOLERANCE 1e-8

// ============
// General Math
// ============

// Constructr
GeneralMath::GeneralMath()
{
    //
}

// Destructor
GeneralMath::~GeneralMath()
{
    //
}

// =============
// Combinatorics
// =============

// Constroctor
Combinatorics::Combinatorics()
{
    //
}

// Destructor
Combinatorics::~Combinatorics()
{
    //
}

// Factorial
unsigned long int Combinatorics::Factorial(unsigned int integer)
{
    if(integer == 0 || integer == 1)
    {
        return 1;
    }

    unsigned long int factorial = 1;
    for(unsigned long int i=2; i<=integer; i++)
    {
        factorial *= i;
    }

    return factorial;
}

// Permutation
unsigned long int Combinatorics::Permutation(unsigned int n, unsigned int k)
{
    if(k>n)
    {
        std::cerr << "k should be less or equal to n." << std::endl;
    }

    if(n == 0 || n == 1 || k == 0)
    {
        return 1;
    }

    unsigned long int permutation = 1;
    for(unsigned long int i=n-k+1; i<=n; i++)
    {
        permutation *= i;
    }
    
    return permutation;
}

// Combination
unsigned long int Combinatorics::Combination(unsigned int n, unsigned int k)
{
    if(k>n)
    {
        std::cerr << "k should be less or equal to n." << std::endl;
    }

    if(n == 0 || n == 1 || n == k || k == 0)
    {
        return 1;
    }

    return Combinatorics::Permutation(n,k) / Combinatorics::Factorial(k);
}

// ======
// Vector
// ======

// ===========
// Constructor
// ===========

Vector::Vector()
{
    //
}

// ==========
// Destructor
// ==========

Vector::~Vector()
{
    //
}

// ==================
// Find Maximum Value
// ==================

double Vector::FindMaximumValue(const double *VectorArray, unsigned int VectorArrayLength)
{
    // Initialize Maxumum Entry
    double MaximumValue = VectorArray[0];

    // Search for larger value
    for(unsigned int VectorArrayIterator = 1; VectorArrayIterator < VectorArrayLength; VectorArrayIterator++)
    {
        // Check with other entries
        if(MaximumValue < VectorArray[VectorArrayIterator])
        {
            // Update maximum value
            MaximumValue = VectorArray[VectorArrayIterator];
        }
    }

    return MaximumValue;
}

// ======
// Matrix
// ======

// ===========
// Constructor
// ===========

Matrix::Matrix()
{
    //
}

// ==========
// Destructor
// ==========

Matrix::~Matrix()
{
    //
}

// =========================
// Determinant Of 2x2 Matrix
// =========================

// Input matrix with 2D array
inline double Matrix::DeterminantOf2x2Matrix(const double Matrix[2][2])
{
    return Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
}

// Input matrix with entries
inline double Matrix::DeterminantOf2x2Matrix(
        const double & Matrix_11, const double & Matrix_12,
        const double & Matrix_21, const double & Matrix_22)
{
    return Matrix_11 * Matrix_22 - Matrix_12 * Matrix_21;
}

// Input matrix with column vectors
inline double Matrix::DeterminantOf2x2Matrix(const double Column1[2],const double Column2[2])
{
    return Column1[0] * Column2[1] - Column2[0] * Column1[1];
}

// =========================
// Determinant Of 3x3 Matrix
// =========================

// Input matrix with 2D array
inline double Matrix::DeterminantOf3x3Matrix(const double Matrix[3][3])
{
    return Matrix[0][0] * (Matrix[1][1] * Matrix[2][2] - Matrix[1][2] * Matrix[2][1]) 
         - Matrix[0][1] * (Matrix[1][0] * Matrix[2][2] - Matrix[1][2] * Matrix[2][0])
         + Matrix[0][2] * (Matrix[1][0] * Matrix[2][1] - Matrix[1][1] * Matrix[2][0]);
}

// Input matrix with entries
inline double Matrix::DeterminantOf3x3Matrix(
        const double & Matrix_11, const double & Matrix_12, const double & Matrix_13,
        const double & Matrix_21, const double & Matrix_22, const double & Matrix_23,
        const double & Matrix_31, const double & Matrix_32, const double & Matrix_33)
{
    return Matrix_11 * (Matrix_22 * Matrix_33 - Matrix_23 * Matrix_32)
         - Matrix_12 * (Matrix_21 * Matrix_33 - Matrix_23 * Matrix_31)
         + Matrix_13 * (Matrix_21 * Matrix_32 - Matrix_22 * Matrix_31);
}

// Input matrix with column vectors
inline double Matrix::DeterminantOf3x3Matrix(
        const double Column1[3], const double Column2[3], const double Column3[3])
{
    return Column1[0] * (Column2[1] * Column3[2] - Column3[1] * Column2[2])
         - Column2[0] * (Column1[1] * Column3[2] - Column3[1] * Column1[2])
         + Column3[0] * (Column1[1] * Column2[2] - Column2[1] * Column1[2]);
}

// =========================
// Determinant Of NxN Matrix
// =========================

// Input matrix with 2D array
inline double Matrix::DeterminantOfNxNMatrix(double **Matrix, unsigned int MatrixSize)
{
    // Convert matrix to row vector //

    // Vector size
    unsigned int VectorizedMatrixSize = MatrixSize * MatrixSize;

    // Declare vector
    double VectorizedMatrix[VectorizedMatrixSize];

    // Set the vector
    for(unsigned int MatrixRowIterator = 0; MatrixRowIterator < MatrixSize; MatrixRowIterator++)
    {
        for(unsigned int MatrixColumnIterator = 0; MatrixColumnIterator < MatrixSize; MatrixColumnIterator++)
        {
            VectorizedMatrix[MatrixRowIterator * MatrixSize + MatrixColumnIterator] = Matrix[MatrixRowIterator][MatrixColumnIterator];
        }
    }

    // Find Determinant
    return Matrix::DeterminantOfNxNMatrix(VectorizedMatrix,MatrixSize);
}

// Input matrix with vectorized array
inline double Matrix::DeterminantOfNxNMatrix(const double *VectorizedMatrix, unsigned int MatrixSize)
{
    // A single number
    if(MatrixSize == 1)
    {
        return *VectorizedMatrix;
    }

    // Declare determinant
    double MatrixDeterminant = 0;

    // Initiliaze Cofactor sign
    int CofactorSign = -1;

    for(unsigned int MatrixColumnIterator = 0; MatrixColumnIterator < MatrixSize; MatrixColumnIterator++)
    {
        // Minor Matrix size
        unsigned int MinorMatrixSize = MatrixSize - 1;

        // Declare Minor Matrix
        double MinorMatrix[MinorMatrixSize][MinorMatrixSize];

        // Construct Minor Matrix //

        // Column terator for minor matrix
        unsigned int MinorColumnIterator = 0;

        // Loop over columns of actual main matrix (again)
        for(unsigned int MatrixColumnIterator2 = 0; MatrixColumnIterator2 < MatrixSize; MatrixColumnIterator2++)
        {
            // Removing target column by skipping
            if(MatrixColumnIterator2 == MatrixColumnIterator)
            {
                continue;
            }

            // Loop over rows of minor matrix
            for(unsigned int MinorRowIterator = 0; MinorRowIterator < MinorMatrixSize; MinorRowIterator++)
            {
                // Equivallent indexing in vectorized array
                unsigned int VectorizedIndex = (MinorRowIterator+1) * MatrixSize + MatrixColumnIterator;

                // Construct minor matrix form main matrix
                MinorMatrix[MinorRowIterator][MinorColumnIterator] = VectorizedMatrix[VectorizedIndex];
            }

            // Update Minor matrix column iterator
            MinorColumnIterator++;
        }

        // Convert Minor matrix to vectorized matrix //

        // Size of vectorized minor matrix
        unsigned int VectorizedMinorSize = MinorMatrixSize * MinorMatrixSize;

        // Declare vectorized minor matrix
        double VectorizedMinorMatrix[VectorizedMinorSize];

        // Convert minor matrix to vectorized form
        for(unsigned int MinorRowIterator = 0; MinorRowIterator < MinorMatrixSize; MinorRowIterator++)
        {
            for(unsigned int MinorColumnIterator = 0; MinorColumnIterator < MinorMatrixSize; MinorColumnIterator++)
            {
                VectorizedMinorMatrix[MinorRowIterator * MinorMatrixSize + MinorColumnIterator] = MinorMatrix[MinorRowIterator][MinorColumnIterator];
            }
        }

        // Determinant of Minor Matrix
        double MinorDeterminant = Matrix::DeterminantOfNxNMatrix(VectorizedMinorMatrix,MinorMatrixSize);

        // Update Cofactor sign
        CofactorSign *= -1;

        // Update Determinant
        MatrixDeterminant += CofactorSign * VectorizedMatrix[MatrixColumnIterator] * MinorDeterminant;
    }

    return MatrixDeterminant;
}

// =================================
// Direct Solve 2x2 Linear Equations
// =================================

void Matrix::DirectSolve2x2LinearEquations(
        const double CoefficientMatrix[2][2],
        const double KnownVector[2],
        double UnknownVector[2])
{
    // Using Cramer's rule //

    // Determinant of Coefficient matrix
    double CoefficientMatrixDeterminant = Matrix::DeterminantOf2x2Matrix(CoefficientMatrix);

    // Determinant of First Cramer's matrix
    double Cramer1Determinant = Matrix::DeterminantOf2x2Matrix(
            KnownVector[0],CoefficientMatrix[0][1],
            KnownVector[1],CoefficientMatrix[1][1]);

    // Solve for first unknown
    UnknownVector[0] = Cramer1Determinant / CoefficientMatrixDeterminant;

    // Determinant of Second Cramer's matrix
    double Cramer2Determinant = Matrix::DeterminantOf2x2Matrix(
            CoefficientMatrix[0][0],KnownVector[0],
            CoefficientMatrix[1][0],KnownVector[1]);

    // Solve for second unknown
    UnknownVector[1] = Cramer2Determinant / CoefficientMatrixDeterminant;
}

// ===========================================
// Direct Solve 2x2 Transpose Linear Equations
// ===========================================

void Matrix::DirectSolve2x2TransposeLinearEquations(
        double TransposedMatrix[2][2],
        double KnownVector[2],
        double UnknownVector[2])
{
    // Columns of the original matrix
    double *Column1 = TransposedMatrix[0];
    double *Column2 = TransposedMatrix[1];

    // Determinant of Transposed matrix
    double TransposedMatrixDeterminant = Matrix::DeterminantOf2x2Matrix(Column1,Column2);

    // Determinant of first Cramer's matrix
    double Cramer1Determinant = Matrix::DeterminantOf2x2Matrix(KnownVector,Column2);

    // Determinant of second Cramer's matrix
    double Cramer2Determinant = Matrix::DeterminantOf2x2Matrix(Column1,KnownVector);

    // Solve for first unknown
    UnknownVector[0] = Cramer1Determinant / TransposedMatrixDeterminant;

    // Solve for second unknown
    UnknownVector[1] = Cramer2Determinant / TransposedMatrixDeterminant;
}

// =================================
// Direct Solve 3x3 Linear Equations
// =================================

void Matrix::DirectSolve3x3LinearEquations(
        const double CoefficientMatrix[3][3],
        const double KnownVector[3],
        double UnknownVector[3])
{
    // Using Cramer's rule //

    // Determinant of coefficient matrix
    double CoefficientMatrixDeterminant = Matrix::DeterminantOf3x3Matrix(CoefficientMatrix);

    // Determinant of first Cramer's matrix
    double Cramer1Determinant = Matrix::DeterminantOf3x3Matrix(
            KnownVector[0],CoefficientMatrix[0][1],CoefficientMatrix[0][2],
            KnownVector[1],CoefficientMatrix[1][1],CoefficientMatrix[1][2],
            KnownVector[2],CoefficientMatrix[2][1],CoefficientMatrix[2][2]);

    // Solve for first unknown
    UnknownVector[0] = Cramer1Determinant / CoefficientMatrixDeterminant;

    // Determinant of second Cramer's matrix
    double Cramer2Determinant = Matrix::DeterminantOf3x3Matrix(
            CoefficientMatrix[0][0],KnownVector[0],CoefficientMatrix[0][2],
            CoefficientMatrix[1][0],KnownVector[1],CoefficientMatrix[1][2],
            CoefficientMatrix[2][0],KnownVector[2],CoefficientMatrix[2][2]);

    // Solve for second unknown
    UnknownVector[1] = Cramer2Determinant / CoefficientMatrixDeterminant;

    // Determinant of third Cramer's matrix
    double Cramer3Determinant = Matrix::DeterminantOf3x3Matrix(
            CoefficientMatrix[0][0],CoefficientMatrix[0][1],KnownVector[0],
            CoefficientMatrix[1][0],CoefficientMatrix[1][1],KnownVector[1],
            CoefficientMatrix[2][0],CoefficientMatrix[2][1],KnownVector[2]);

    // Solve for third unknown
    UnknownVector[2] = Cramer3Determinant / CoefficientMatrixDeterminant;
}

// ===========================================
// Direct Solve 3x3 Transpose Linear Equations
// ===========================================

void Matrix::DirectSolve3x3TransposeLinearEquations(
        double TransposedMatrix[3][3],
        double KnownVector[3],
        double UnknownVector[3])
{
    // Columns of the original matrix
    double *Column1 = TransposedMatrix[0];
    double *Column2 = TransposedMatrix[1];
    double *Column3 = TransposedMatrix[2];

    // Determinant of Transposed matrix
    double TransposedMatrixDeterminant = Matrix::DeterminantOf3x3Matrix(Column1,Column2,Column3);

    // Determinant of first Cramer's matrix
    double Cramer1Determinant = Matrix::DeterminantOf3x3Matrix(KnownVector,Column2,Column3);

    // Determinant of second Cramer's matrix
    double Cramer2Determinant = Matrix::DeterminantOf3x3Matrix(Column1,KnownVector,Column3);

    // Determinant of third Cramer's matrix
    double Cramer3Determinant = Matrix::DeterminantOf3x3Matrix(Column1,Column2,KnownVector);

    // Solve for first unknown
    UnknownVector[0] = Cramer1Determinant / TransposedMatrixDeterminant;

    // Solve for second unknown
    UnknownVector[1] = Cramer2Determinant / TransposedMatrixDeterminant;

    // Solve for third unknown
    UnknownVector[2] = Cramer3Determinant / TransposedMatrixDeterminant;
}

// =================================
// Direct Solve NxN Linear Equations
// =================================

void Matrix::DirectSolveNxNLinearEquations(
        double **CoefficientMatrix,
        double *KnownVector,
        unsigned int NumberOfEquations,
        double *UnknownVector)
{
    // Using Crame's rule

    // Determinant of coefficient matrix
    double CoefficientMatrixDeterminant = Matrix::DeterminantOfNxNMatrix(CoefficientMatrix,NumberOfEquations);

    // Check small determinants
    if(CoefficientMatrixDeterminant < TOLERANCE)
    {
        std::cerr << "Coefficient matrix is singular." << std::endl;
    }

    // Find transpose of coefficient matrix
    double TransposeMatrix[NumberOfEquations][NumberOfEquations];

    for(unsigned int RowIterator = 0; RowIterator < NumberOfEquations; RowIterator++)
    {
        for(unsigned int ColumnIterator = 0; ColumnIterator < NumberOfEquations; ColumnIterator++)
        {
            TransposeMatrix[RowIterator][ColumnIterator] = CoefficientMatrix[ColumnIterator][RowIterator];
        }
    }

    // Apply Cramer's rule for each column of matrix
    for(unsigned int ColumnIterator = 0; ColumnIterator < NumberOfEquations; ColumnIterator++)
    {
        // Declare Cramer's Matrix
        double *CramerMatrix[NumberOfEquations];

        // Set the Crame'r Matrix'
        for(unsigned int ColumnIterator2 = 0; ColumnIterator < NumberOfEquations; ColumnIterator++)
        {
            if(ColumnIterator2 == ColumnIterator)
            {
                CramerMatrix[ColumnIterator2] = KnownVector;
            }
            else
            {
                // Get thecolumn of coefficient matrix by picking row of transpose matrix
                CramerMatrix[ColumnIterator2] = TransposeMatrix[ColumnIterator];
            }
        }

        // Find determinant of Cramer matrix
        double CramerMatrixDeterminant = Matrix::DeterminantOfNxNMatrix(CramerMatrix,NumberOfEquations);

        // Solve the unknown value by cramer's rule
        UnknownVector[ColumnIterator] = CramerMatrixDeterminant / CoefficientMatrixDeterminant;
    }
}

// ===========================================
// Direct Solve NxN Transpose Linear Equations
// ===========================================

void Matrix::DirectSolveNxNTransposeLinearEquations(
        double **TransposedMatrix,
        double *KnownVector,
        unsigned int NumberOfEquations,
        double *UnknownVector)
{

}

// ====================================
// Iterative Solve NxN Linear Equations
// ====================================

void Matrix::IterativeSolveNxNLinearEquations(
        const double **CoefficientMatrix,
        const double *KnownVector,
        unsigned int NumberOfEquations,
        double *UnknownVector)
{
    // Gauss Seidel Method //

    // Initialize Maximum Relative Error
    double MaximumRelativeError = 1+TOLERANCE;

    // Maximum number Of iterations
    unsigned int MaximumNumberOfIterations = 100;

    // Counter number of Iterations
    unsigned int NumberOfIterations = 0;

    // Declare Previous unkown vector
    double PreviousUnknownVector[NumberOfEquations];

    // Initialize Previous unknown vector as a guess
    for(unsigned int VectorIndex = 0; VectorIndex < NumberOfEquations; VectorIndex++)
    {
        PreviousUnknownVector[VectorIndex] = 0;
        UnknownVector[VectorIndex] = 0;
    }

    // Convergence iteration
    while(MaximumRelativeError < TOLERANCE && NumberOfIterations < MaximumNumberOfIterations)
    {

        // Solve for each entry
        for(unsigned int i = 0; i < NumberOfEquations; i++)
        {
            // Initialize updated vector
            UnknownVector[i] = KnownVector[i];

            // Multiple subtraction
            for(unsigned int j = 0; j < NumberOfEquations; j++)
            {
                if(i == j)
                {
                    continue;
                }
                else
                {
                    UnknownVector[i] -= CoefficientMatrix[i][j] * UnknownVector[j];
                }
            }

            // Divide by diagonal element
            if(CoefficientMatrix[i][i] != 0)
            {
                UnknownVector[i] /= CoefficientMatrix[i][i];
            }
            else
            {
                std::cerr << "Division by zero. Matrix needs pivoting." << std::endl;
            }
        }

        // Absolute error
        double AbsoluteError[NumberOfEquations];

        for(unsigned int VectorIndex = 0; VectorIndex < NumberOfEquations; VectorIndex++)
        {
            AbsoluteError[VectorIndex] = fabs(UnknownVector[VectorIndex] - PreviousUnknownVector[VectorIndex]);
        }

        // Find maximum value
        double MaximumOfUnknownVector = Vector::FindMaximumValue(PreviousUnknownVector,NumberOfEquations);

        // Find maximum of absolute Error
        double MaximumAbsoluteError = Vector::FindMaximumValue(AbsoluteError,NumberOfEquations);

        // Update maximum Relative Error
        MaximumRelativeError = MaximumAbsoluteError / (MaximumOfUnknownVector + TOLERANCE);

        // Update Number of iteratios
        NumberOfIterations++;

        // Update Previous Unknown vectors
        for(unsigned int VectorIndex = 0; VectorIndex < NumberOfEquations; VectorIndex++)
        {
            PreviousUnknownVector[VectorIndex] = UnknownVector[VectorIndex];
        }
    }
}

// ==========
// Polynomial
// ==========

// ============
// Constructors
// ============

Polynomial::Polynomial()
{
    this->Coefficients = new RationalNumber[1];
    this->Coefficients[0] = 0;
    this->Degree = 0;
}

Polynomial::Polynomial(int integer)
{
    this->Coefficients = new RationalNumber[1];
    this->Coefficients[0] = integer;
    this->Degree = 0;
}

Polynomial::Polynomial(double Double)
{
    this->Coefficients = new RationalNumber[1];
    this->Coefficients[0] = Double;
    this->Degree = 0;
}

Polynomial::Polynomial(RationalNumber rational)
{
    this->Coefficients = new RationalNumber[1];
    this->Coefficients[0] = rational;
    this->Degree = 0;
}

Polynomial::Polynomial(const Polynomial & rhs)
{
    // Initialize (it is necessary here)
    this->Coefficients = new RationalNumber[1];
    this->Coefficients[0] = 0;
    this->Degree = 0;

    // Replace with new values from rhs
    this->SetCoefficients(rhs.GetCoefficients(),rhs.GetDegree());
    this->Degree = rhs.GetDegree();
    this->Squeeze();
}

Polynomial::~Polynomial()
{
    delete [] this->Coefficients;
    this->Coefficients = NULL;
}

// ===================
// Accessors, Mutators
// ===================

// Get Coefficients
RationalNumber * Polynomial::GetCoefficients() const
{
    return this->Coefficients;
}

// Get Coefficient
RationalNumber Polynomial::GetCoefficient(unsigned int index) const
{
    return this->Coefficients[index];
}

// Set Coefficient
void Polynomial::SetCoefficient(unsigned int index, double coef)
{
    this->SetCoefficient(index,RationalNumber(coef));
}

// Set Coefficient
void Polynomial::SetCoefficient(unsigned int index,RationalNumber RNcoef)
{
    if(index<=this->Degree)
    {
        this->Coefficients[index] = RNcoef;
    }
    else
    {
        this->Extend(index);
        this->Coefficients[index] = RNcoef;
    }
}

// Set Coefficients
void Polynomial::SetCoefficients(double * DoubleCoefs, unsigned int degree)
{
    // Convert integers to rational number
    RationalNumber *RNcoefs = new RationalNumber[degree+1];
    for(unsigned int i=0; i<=degree; i++)
    {
        RNcoefs[i] = RationalNumber(DoubleCoefs[i]);
    }

    // Pass to overloaded method with rational number
    this->SetCoefficients(RNcoefs,degree);

    // freed memory
    delete [] RNcoefs;
    RNcoefs = NULL;
}

// Set Coefficients
void Polynomial::SetCoefficients(RationalNumber * RNCoefs, unsigned int degree)
{
    // Remove old allocation
    if(this->Coefficients != NULL)
    {
        delete [] this->Coefficients;
        this->Coefficients = NULL;
    }

    // Allocate new coefficients
    this->Degree = degree;
    this->Coefficients = new RationalNumber[degree+1];
    
    // Copy coefficients
    for(unsigned int i=0; i<=degree; i++)
    {
        this->Coefficients[i] = RNCoefs[i];
    }
}

// Add Next Coefficient
void Polynomial::AddNextCoefficient(double coef)
{
    this->AddNextCoefficient(RationalNumber(coef));

}

// Add Next Coefficient
void Polynomial::AddNextCoefficient(RationalNumber RNcoef)
{
    unsigned short int Degree_old = this->GetDegree();
    this->Extend(Degree_old +1);
    this->Coefficients[Degree_old+1] = RNcoef;
}

// Get Degree
unsigned int Polynomial::GetDegree() const
{
    return this->Degree;
}

// Set Degree
void Polynomial::SetDegree(unsigned int degree)
{
    this->Resize(degree);
}

// Get As String
const char * Polynomial::GetAsString() const
{
    std::stringstream ss;
    for(unsigned int i=0; i<=this->GetDegree(); i++)
    {
        RationalNumber r = this->GetCoefficient(i);
        ss << std::string(r.GetAsString());
        if(i < this->GetDegree())
        {
            ss << ", ";
        }
    }
    return ss.str().c_str();
}

// =======
// Methods
// =======

// Differentiate
Polynomial Polynomial::Differentiate()
{
    return Polynomial::Differentiate(*this,1);
}

// Static Differentiate
Polynomial Polynomial::Differentiate(const Polynomial & P)
{
    return Polynomial::Differentiate(P,1);
}

// Differentiate n-th
Polynomial Polynomial::Differentiate(unsigned int order)
{
    return Polynomial::Differentiate(*this,order);
}

// Static Differentiate n-th
Polynomial Polynomial::Differentiate(const Polynomial & P, unsigned int order)
{
    if(order == 0)
    {
        return P;
    }

    // Differentiate of P
    Polynomial dP;
    dP.SetDegree(P.GetDegree()-order);
    for(unsigned int i=0; i<=dP.GetDegree(); i++)
    {
        RationalNumber coef = P.GetCoefficient(i+order) * Combinatorics::Permutation(i+order,order);
        dP.SetCoefficient(i,coef);
    }

    dP.Squeeze();

    return dP;
}

// Integrate
Polynomial Polynomial::Integrate()
{
    return Polynomial::Integrate(*this,1);
}

// Static Integrate
Polynomial Polynomial::Integrate(const Polynomial & P)
{
    return Polynomial::Integrate(P,1);
}

// Integrate n-th
Polynomial Polynomial::Integrate(unsigned int order)
{
    return Polynomial::Integrate(*this,order);
}

// Static Integrate n-th
Polynomial Polynomial::Integrate(const Polynomial & P, unsigned int order)
{
    if(order == 0)
    {
        return P;
    }

    // Integration of P
    Polynomial IntP;
    IntP.SetDegree(P.GetDegree()+order);
    for(unsigned int i=order; i<=IntP.GetDegree(); i++)
    {
        // Coeffcient generated due to integration
        RationalNumber IntegrationCoef(1);
        for(unsigned int j=1; j<=order; j++)
        {
            IntegrationCoef *= RationalNumber(1,i-order+j);
        }

        // multiplying to corfficient of P
        RationalNumber coef = P.GetCoefficient(i-order) * IntegrationCoef;

        // Set coefficient of IntP
        IntP.SetCoefficient(i,coef);
    }

    IntP.Squeeze();

    return IntP;
}

// Evaluate
RationalNumber Polynomial::Evaluate(const RationalNumber & rational)
{
    return Polynomial::Evaluate(*this,rational);
}

// Static Evaluate
RationalNumber Polynomial::Evaluate(const Polynomial & P, const RationalNumber & rational)
{
    RationalNumber value(0);
    for(int i=0; i<=P.GetDegree(); i++)
    {
        value += P.GetCoefficient(i) * (rational^i);
    }
    return value.NormalForm();
}

// Resize
void Polynomial::Resize(unsigned int NewDegree)
{
    if(NewDegree < this->GetDegree())
    {
        this->Squeeze(NewDegree);
    }
    else if(NewDegree > this->GetDegree())
    {
        this->Extend(NewDegree);
    }
}

// Squeeze
void Polynomial::Squeeze()
{
    // Find zero coefficient with smallest index
    unsigned int VanishedIndex;
    bool FoundZero = false;
    for(unsigned int i=this->GetDegree(); i>=0; i--)
    {
        if(this->Coefficients[i] == 0)
        {
            FoundZero = true;
            VanishedIndex = i;
        }
        else
        {
            break;
        }
    }

    // Squeeze
    if(FoundZero == true && VanishedIndex > 0)
    {
        this->Squeeze(VanishedIndex-1);
    }
}

// Squeeze
void Polynomial::Squeeze(unsigned int NewDegree)
{
    // Check if new degree is less
    if(NewDegree > this->GetDegree())
    {
        std::cerr << "ERROR: Can not squeeze to greater degree." << std::endl;
    }

    // Squeeze
    if(NewDegree < this->GetDegree())
    {
        // Save previous data
        unsigned int Degree_old = this->GetDegree();
        RationalNumber *Coefficients_old = new RationalNumber[NewDegree+1];
        for(unsigned int i=0; i<=NewDegree; i++)
        {
            Coefficients_old[i] = this->Coefficients[i];
        }

        // Update coefficients array
        this->Degree = NewDegree;
        delete [] this->Coefficients;
        this->Coefficients = new RationalNumber[NewDegree+1];
        for(unsigned int i=0; i<=NewDegree; i++)
        {
            this->Coefficients[i] = Coefficients_old[i];
        }

        // freed memory
        delete [] Coefficients_old;
        Coefficients_old = NULL;
    }
}

// Extend
void Polynomial::Extend(unsigned int NewDegree)
{
    // Check if new degree if greater
    if(NewDegree < this->GetDegree())
    {
        std::cerr << "ERROR: Can not extend polynomial to smaller degree." << std::endl;
    }

    // Extend
    if(NewDegree > this->GetDegree())
    {
        // Store old coefficients
        unsigned int Degree_old = this->GetDegree();
        RationalNumber *Coefficients_old = new RationalNumber[Degree_old+1];
        for(unsigned int i=0; i<=Degree_old; i++)
        {
            Coefficients_old[i] = this->Coefficients[i];
        }

        // Extend member data array
        delete [] this->Coefficients;
        this->Coefficients = NULL;
        this->Coefficients = new RationalNumber[NewDegree+1];
        this->Degree = NewDegree;

        // Retrive old data
        for(unsigned int i=0; i<=Degree_old; i++)
        {
            this->Coefficients[i] = Coefficients_old[i];
        }

        // Fill rest of extention with zeros
        for(unsigned int i=Degree_old+1; i<=NewDegree; i++)
        {
            this->Coefficients[i] = 0;
        }

        // freed memory
        delete [] Coefficients_old;
        Coefficients_old = NULL;
    }
}

// Print
void Polynomial::Print() const
{
    Polynomial::Print(std::cout,*this);
}

// Static Print
void Polynomial::Print(std::ostream & os, const Polynomial & P)
{
    unsigned int degree = P.GetDegree();
    for(unsigned int i=0; i<=degree; i++)
    {
        os << P.GetCoefficient(i).GetAsString();
        if(i < degree)
        {
            os << ", ";
        }
    }
    os << std::endl;
}

// =========
// Operators
// =========

// +
Polynomial Polynomial::operator+(const Polynomial & rhs) const
{
    // Determine min/max degree polynomials
    unsigned int ThisDegree = this->GetDegree();
    unsigned int rhsDegree = rhs.GetDegree();
    unsigned int MinDegree = (ThisDegree > rhsDegree ? rhsDegree : ThisDegree);
    unsigned int MaxDegree = (ThisDegree > rhsDegree ? ThisDegree : rhsDegree);

    // Sum polynomial
    Polynomial P;
    P.SetDegree(MaxDegree);

    // Adding coefficients of common degrees
    for(unsigned int i=0; i<= MinDegree; i++)
    {
        P.SetCoefficient(i,this->GetCoefficient(i) + rhs.GetCoefficient(i));
    }

    // Adding coefficients of higher degrees
    if(ThisDegree > rhsDegree)
    {
        for(unsigned int i=MinDegree+1; i <= MaxDegree; i++)
        {
            P.SetCoefficient(i,this->GetCoefficient(i));
        }
    }
    else if(ThisDegree < rhsDegree)
    {
        for(unsigned int i=MinDegree+1; i <= MaxDegree; i++)
        {
            P.SetCoefficient(i,rhs.GetCoefficient(i));
        }
    }

    P.Squeeze();

    return P;
}

// -
Polynomial Polynomial::operator-(const Polynomial & rhs) const
{
    return *this + (rhs * (-1));
}

// *
Polynomial Polynomial::operator*(const Polynomial & rhs) const
{
    // Left and right of product
    unsigned int p = this->GetDegree();
    unsigned int q = rhs.GetDegree();

    // Product Polynomial
    Polynomial Product;
    Product.SetDegree(p + q);

    // Convolution
    for(unsigned int i=0; i<= Product.GetDegree(); i++)
    {
        // Index bounds
        unsigned int j_min = (0 < static_cast<signed int>(i-p) ? static_cast<signed int>(i-p) : 0); // max of 0 or i-q
        unsigned int j_max = (i < q ? i : q); // min of q or i

        RationalNumber r(0);
        for(unsigned int j=j_min; j<=j_max; j++)
        {
            r += rhs.GetCoefficient(j) * this->GetCoefficient(i-j);
        }

        Product.SetCoefficient(i,r);
    }

    Product.Squeeze();

    return Product;
}

// ^
Polynomial Polynomial::operator^(const int power) const
{
    if(power == 1)
    {
        return *this;
    }

    Polynomial P(*this);
    for(unsigned int i=2; i<=power; i++)
    {
        P *= *this;
    }

    return P;
}

// =
Polynomial & Polynomial::operator=(const Polynomial & rhs)
{
    if(this == &rhs)
    {
        return *this;
    }

    // Update this
    this->SetCoefficients(rhs.GetCoefficients(),rhs.GetDegree());

    return *this;
}

// +=
Polynomial & Polynomial::operator+=(const Polynomial & rhs)
{
    // Addition
    Polynomial P = *this + rhs;

    // Update this
    this->SetCoefficients(P.GetCoefficients(),P.GetDegree());

    return *this;
}

// -=
Polynomial & Polynomial::operator-=(const Polynomial & rhs)
{
    // Subtraction
    Polynomial P = *this - rhs;

    // Update this
    this->SetCoefficients(P.GetCoefficients(),P.GetDegree());

    return *this;
}

// *=
Polynomial & Polynomial::operator*=(const Polynomial & rhs)
{
    // Multiplication
    Polynomial P = *this * rhs;

    // Update this
    this->SetCoefficients(P.GetCoefficients(),P.GetDegree());

    return *this;
}

// []
RationalNumber Polynomial::operator[](unsigned int offset) const
{
    return this->GetCoefficient(offset);
}

// ==
bool Polynomial::operator==(const Polynomial & rhs) const
{
    if(this->GetDegree() != rhs.GetDegree())
    {
        return false;
    }

    for(unsigned int i=0; i<=this->GetDegree(); i++)
    {
        if(this->GetCoefficient(i) != rhs.GetCoefficient(i))
        {
            return false;
        }
    }

    return true;
}

// !=
bool Polynomial::operator!=(const Polynomial & rhs) const
{
    return !(this->operator==(rhs));
}

// <<
std::ostream & operator<<(std::ostream & os, const Polynomial & rhs)
{
    unsigned int degree = rhs.GetDegree();
    for(unsigned int i=0; i<=degree; i++)
    {
        os << rhs.GetCoefficient(i).GetAsString();
        if(i < degree)
        {
            os << ", ";
        }
    }
    return os;
}

// Converstion to int
Polynomial::operator int() const
{
    return static_cast<int>(this->GetCoefficient(0).IntegerPart());
}

// Conversion to double
Polynomial::operator double() const
{
    return this->GetCoefficient(0).Evaluate();
}

// Converstion to Rational Number
Polynomial::operator RationalNumber() const
{
    return this->GetCoefficient(0);
}

// ===============
// Rational Number
// ===============

// ============
// Constructors
// ============

// Default
RationalNumber::RationalNumber():
    Numerator(0),
    Denominator(1)
{
}

// With double or integer
RationalNumber::RationalNumber(double Double)
{
    double PositiveDouble = fabs(Double);
    if(fabs(PositiveDouble - static_cast<int>(PositiveDouble)) < 1e-15)
    {
        this->Numerator = static_cast<int>(Double);
        this->Denominator = 1;
        this->AdjustSign();
    }
    else
    {
        // Find number of fractional part digits
        std::stringstream ss;
        ss << fabs(PositiveDouble-static_cast<int>(PositiveDouble));
        unsigned short int FractionalDigits = ss.str().length() - 2;
        unsigned short int MaxPrecision = 8;
        unsigned short int DigitsPrecision = (FractionalDigits < MaxPrecision ? FractionalDigits : MaxPrecision);

        // Convert double to rational number
        this->Numerator = round(Double * pow(10,DigitsPrecision));
        this->Denominator = pow(10,DigitsPrecision);
        this->AdjustSign();
        this->Simplify();
    }
}

// Full definition
RationalNumber::RationalNumber(long int numerator, long int denominator):
    Numerator(numerator),
    Denominator(denominator)
{
    this->CheckValid();
    this->AdjustSign();
    this->Simplify();
}

// Copy Constructor
RationalNumber::RationalNumber(const RationalNumber & rhs)
{
    this->Numerator = rhs.GetNumerator();
    this->Denominator = rhs.GetDenominator();
    this->AdjustSign();
    this->Simplify();
}

// Destructor
RationalNumber::~RationalNumber()
{
}

// ===================
// Accessors, Mutators
// ==================

// Set Denominator
void RationalNumber::SetDenominator(long int denominator)
{
    this->Denominator = denominator;
    this->CheckValid();
}

// Get As String
const char * RationalNumber::GetAsString() const
{
    std::stringstream ss;
    ss << "(" << this->Numerator << "/" << this->Denominator << ")";
    return ss.str().c_str();
}

// ==============
// Public Methods
// ==============

// Integer Part
long int RationalNumber::IntegerPart() const
{
    return (this->Numerator - (this->Numerator % this->Denominator)) / this->Denominator;
}

// Rational Part
RationalNumber RationalNumber::FractionalPart() const
{
    long int numerator = this->Numerator % this->Denominator;
    return RationalNumber(numerator,this->Denominator);
}

// Inverse
RationalNumber RationalNumber::Inverse() const
{
    return RationalNumber(this->Denominator,this->Numerator);
}

// Opposite
RationalNumber RationalNumber::Opposite() const
{
    return RationalNumber(-this->Numerator,this->Denominator);
}

// Greatest Common Divisor
long int RationalNumber::GreatestCommonDivisor(long int a, long int b)
{
    if(a<0)
    {
        a = -a;
    }

    if(b<0)
    {
        b = -b;
    }

    for(;;)
    {
        if(a == 0) return b;
        b %= a;
        if(b == 0) return a;
        a %= b;
    }
}

// Least Common Multiple
long int RationalNumber::LeastCommonMultiple(long int a, long int b)
{
    long int gcd = RationalNumber::GreatestCommonDivisor(a,b);
    return gcd ? ((a/gcd)*b) : 0;
}

// Simplify
RationalNumber & RationalNumber::Simplify()
{
    long int gcd = RationalNumber::GreatestCommonDivisor(this->Numerator,this->Denominator);
    this->Numerator /= fabs(gcd);
    this->Denominator /= fabs(gcd);
    return *this;
}

// Normal Form
RationalNumber & RationalNumber::NormalForm()
{
    this->CheckValid();
    this->AdjustSign();
    this->Simplify();
    return *this;
}

// Evaluate
double RationalNumber::Evaluate() const
{
    return static_cast<double>(this->Numerator) / static_cast<double>(this->Denominator);
}

// Print
void RationalNumber::Print() const
{
    RationalNumber::Print(std::cout,*this);
}

// Static Print
void RationalNumber::Print(std::ostream & os,const RationalNumber & RN)
{
    os << RN.GetAsString() << std::endl;
}

// =========
// Operators
// =========

// +
RationalNumber RationalNumber::operator+(const RationalNumber & rhs) const
{
    long int lcm = RationalNumber::LeastCommonMultiple(this->Denominator,rhs.GetDenominator());
    long int numerator = (this->Numerator * (lcm/this->Denominator)) + (rhs.GetNumerator() * (lcm/rhs.GetDenominator()));
    RationalNumber tempRN(numerator,lcm);
    return tempRN.Simplify();
}

// -
RationalNumber RationalNumber::operator-(const RationalNumber & rhs) const
{
    return *this + rhs.Opposite();
}

// *
RationalNumber RationalNumber::operator*(const RationalNumber & rhs) const
{
    RationalNumber tempRN(
            this->Numerator * rhs.GetNumerator(),
            this->Denominator * rhs.GetDenominator());
    return tempRN.Simplify();
}

// /
RationalNumber RationalNumber::operator/(const RationalNumber & rhs) const
{
    return *this * rhs.Inverse();
}

// ^
RationalNumber RationalNumber::operator^(int power) const
{
    if(power < 0)
    {
        return RationalNumber(pow(this->Denominator,fabs(power)),pow(this->Numerator,fabs(power)));
    }
    else if(power == 0)
    {
         return RationalNumber(1);
    }
    else
    {
        return RationalNumber(pow(this->Numerator,power),pow(this->Denominator,power));
    }
}

// =
RationalNumber & RationalNumber::operator=(const RationalNumber & rhs)
{
    if(this == &rhs)
    {
        return *this;
    }

    this->Numerator = rhs.GetNumerator();
    this->Denominator = rhs.GetDenominator();

    return *this;
}

// +=
RationalNumber & RationalNumber::operator+=(const RationalNumber & rhs)
{
    RationalNumber tempRN = *this + rhs;
    tempRN.Simplify();
    this->Numerator = tempRN.GetNumerator();
    this->Denominator = tempRN.GetDenominator();
    return *this;
}

// -=
RationalNumber & RationalNumber::operator-=(const RationalNumber & rhs)
{
    return *this += rhs.Opposite();
}

// *=
RationalNumber & RationalNumber::operator*=(const RationalNumber & rhs)
{
    RationalNumber tempRN = *this * rhs;
    tempRN.Simplify();
    this->Numerator = tempRN.GetNumerator();
    this->Denominator = tempRN.GetDenominator();
    return *this;
}

// /=
RationalNumber & RationalNumber::operator/=(const RationalNumber & rhs)
{
    return *this *= rhs.Inverse();
}

// <
bool RationalNumber::operator<(const RationalNumber & rhs) const
{
    if(this->Numerator * rhs.GetDenominator() < this->Denominator * rhs.GetNumerator())
    {
        return true;
    }
    return false;
}

// >
bool RationalNumber::operator>(const RationalNumber & rhs) const
{
    if(this->Numerator * rhs.GetDenominator() > this->Denominator * rhs.GetNumerator())
    {
        return true;
    }
    return false;
}

// <=
bool RationalNumber::operator<=(const RationalNumber & rhs) const
{
    if(this->Numerator * rhs.GetDenominator() <= this->Denominator * rhs.GetNumerator())
    {
        return true;
    }
    return false;
}

// >=
bool RationalNumber::operator>=(const RationalNumber & rhs) const
{
    if(this->Numerator * rhs.GetDenominator() >= this->Denominator * rhs.GetNumerator())
    {
        return true;
    }
    return false;
}

// ==
bool RationalNumber::operator==(const RationalNumber & rhs) const
{
    if(this->GetNumerator() * rhs.GetDenominator() == this->Denominator * rhs.GetNumerator())
    {
        return true;
    }
    return false;
}

// !=
bool RationalNumber::operator!=(const RationalNumber & rhs) const
{
    if(this->Numerator * rhs.GetDenominator() != this->Denominator * rhs.GetNumerator())
    {
        return true;
    }
    return false;
}

// <<
std::ostream & operator<<(std::ostream & os, const RationalNumber & rhs)
{
    os << rhs.GetAsString();
    return os;
}

// Conversion to int
RationalNumber::operator int() const
{
    return this->IntegerPart();
}

// Conversion to double
RationalNumber::operator double() const
{
    return this->Evaluate();
}

// =================
// Protected Methods
// =================

// Check Valid
bool RationalNumber::CheckValid() const
{
    if(this->Denominator == 0)
    {
        std::cerr << "Zero denominator is not allowed." << std::endl;
        return false;
    }
    return true;
}

// Adjust Sign
void RationalNumber::AdjustSign()
{
    if(this->Denominator < 0)
    {
        this->Denominator = -this->Denominator;
        this->Numerator = -this->Numerator;
    }
}
