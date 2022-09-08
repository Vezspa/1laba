#include "Decomposition.h"
#include "Matrix.h"
#include <cassert>
#include "Solver.h"


// 1) Constructors:
Decomposition::Decomposition() {
    size = 0;
    LU = Matrix();
}

Decomposition::Decomposition(const Matrix& any)
{
    // 0. Checking of sizes. If matrix isn't square then error out!
    assert((any.get_cSize() == any.get_rSize()) && "ERROR_MATRIX_IS_NOT_SQUARE");
    //assert((any.get_cSize() == 0) && "ERROR_MATRIX_IS_EMPTY");

    // 1. The data is set there:
    this->size = any.get_cSize();

    LU = Matrix(size, size);

    for (size_t i = 0; i < size; i++)
    {
        LU.at(0, i) = any.at(0, i);
    }
    for (size_t i = 1; i < size; i++)
    {
        LU.at(i, 0) = any.at(i, 0)/any.at(0,0); 
    }
    for (size_t i = 1; i < size; i++)
    {
        for (size_t k = i; k < size; k++)
        {
            double res = 0;
            for (size_t p = 0; p < i-1; p++)
            {
                res += LU.at(i, p) * LU.at(p, k);
            }
            LU.at(i, k) = any.at(i, k) - res;
        }
        for (size_t k = i+1; k < size; k++)
        {
            double res = 0;
            for (size_t p = 0; p < i - 1; p++)
            {
                res += LU.at(k, p) * LU.at(p, i);
            }
            LU.at(k, i) = (any.at(k, i) - res) / LU.at(i, i);
        }
    }
    // По умолчанию: (L\U) = { 10.0, -0.1, 1.0, 10.1 } - LU разложение для матрицы А = {10.0, -1.0, 1.0, 10.0}
}

// 2) Destructor:
Decomposition::~Decomposition()
{
    this->values.clear();
    this->values.shrink_to_fit();
}


// 3) Geters and seters:
const double Decomposition::get_elemL(unsigned int row, unsigned int col) const
{
    // 0. Checking of the indexes!
    assert(((row < this->size) && (col < this->size)) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");

    Matrix L = get_L();

    return L.at(row, col);
}

const double Decomposition::get_elemU(unsigned int row, unsigned int col) const
{
    // 0. Checking of the indexes!
    assert(((row < this->size) && (col < this->size)) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");

    Matrix U = get_U();

    return U.at(row, col);
}

const double Decomposition::get_size() const
{
    return this->size;
}

const Matrix Decomposition::get_L() const
{
    Matrix L(size, size);
    for (size_t i = 0; i < size; i++)
    {
        for (size_t k = 0; k < size; k++)
        {
            if (i > k) { L.at(i, k) = LU.at(i, k); continue; }
            if (i == k) { L.at(i, k) = 1; continue; }
            if (i < k) { L.at(i, k) = 0; }
        }
    }
    return L;
}

const Matrix Decomposition::get_U() const
{
    Matrix U(size, size);
    for (size_t i = 0; i < size; i++)
    {
        for (size_t k = 0; k < size; k++)
        {
            if (i <= k) { U.at(i, k) = LU.at(i, k); continue; }
            if (i > k) { U.at(i, k) = 0; }
        }
    }
    return U;
}