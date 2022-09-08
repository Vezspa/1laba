#include "Matrix.h"
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>

// Stuff functions:
unsigned int get_row2swap(const unsigned int index_diag, const Matrix& Any)
// ������� ������������ ����� ������ ������ �� ������� ���������� �������� ���������������,
// � ������� ������������ ��������� (����� ���������� �������� � ������� � ������� index_diag).
// 
// ���� ����������� �������� �� ������� index2swap ��������� ���������� ���� � ������� =>
// => ������������ ������� ����!
{
	unsigned int index2swap = index_diag;
	while ((index2swap < Any.get_rSize()) && (Any.at(index2swap, index_diag) == 0.0))
	{			
		index2swap = index2swap + 1;
	}

	return index2swap;
};

bool swap_rows(const unsigned int index_diag, Matrix& Any)
// ������� ���������� ���� bool:
//		true, ���� ���� ��������� ������������ �����
//		false, ���� ������������ ����� �� ����
{
	bool swap_flag = false;

	const unsigned int index2swap = get_row2swap(index_diag, Any);

	if (index2swap != Any.get_rSize())
	{
		// ������������ �����:
		swap_flag = true;

		// ������������ �����:
		double buffer_value;
		for (size_t col = index_diag; col < Any.get_cSize(); col++)
		{
			buffer_value = Any.at(index_diag, col);
			Any.at(index_diag, col) = Any.at(index2swap, col);
			Any.at(index2swap, col) = buffer_value;
		}
	}

	return swap_flag;
}

void column_reset(const unsigned int index_diag, Matrix& Any) 
{
	double swaped_value;
	for (size_t row = index_diag + 1; row < Any.get_rSize(); row++)
	{
		swaped_value = Any.at(row, index_diag);
		Any.at(row, index_diag) = 0.0;

		for (size_t col = index_diag + 1; col < Any.get_cSize(); col++) // ���� �������� � ��������. ��������� ������ ��� col = row, �.�. ������ ��� ���������� �������� �� + 1
		{
			//add(Any, row, col, -swaped_value * Any.get_elem(index_diag, col));
			Any.at(row, col) -= swaped_value * Any.at(index_diag, col);
		}
	}
}

// ���������� �������, ������� ����� ������ �� ������������ ������� � ������������ �������� �������������
void row_sub(const unsigned int index_diag, double& det_value, Matrix& Copy)
{
	const double value_diag = Copy.at(index_diag, index_diag);
	if (value_diag != 1.0)
	{
		det_value = det_value * value_diag;

		// ������� ������ �� ������������ �������:
		for (size_t col = index_diag; col < Copy.get_cSize(); col++)
		{
			Copy.at(index_diag, col) /= value_diag;
		}
	}
}


// -1) The private geter gets a linear index:
unsigned int Matrix::get_index(unsigned int row, unsigned int col) const
{
	// n = i - 1 + (j - 1) * rown    � ������, ���� i in [1, rown], � j in [1, coln]  =>
	// => return row - 1 + (col - 1) * this->rown;
	// n = i + j * rown		� ������, ���� i in [0, rown-1], � j in [0, coln-1] =>
	// => return row + col * this->rown;
	// ������ ��������� �������� � exel ����� �������.
	
	//assert((col < this->coln) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE"); // assert(bool = true)
	//assert((row < this->rown) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");

	return row + col * this->rown;
}

// 1) �onstructors:
Matrix::Matrix() : values({10.0, -1.0, 1.0, 10.0}), rown(2), coln(2) {}

Matrix::Matrix(unsigned int rown, unsigned int coln) : coln(coln), rown(rown), values(coln * rown) {}

// 2) Destructior:
Matrix::~Matrix()
{
	values.clear();
	values.shrink_to_fit();
}

// 3) Geters and seters:
const unsigned int Matrix::get_rSize() const
{
	return this->rown;
}

const unsigned int Matrix::get_cSize() const
{
	return this->coln;
}

double& Matrix::at(unsigned int row, unsigned int col)
{
	assert((row < this->get_rSize()) && (col < this->get_cSize()) && "KEK");
	return this->values.at(get_index(row, col));
}

const double& Matrix::at(unsigned int row, unsigned int col) const
{
	assert((row < this->get_rSize()) && (col < this->get_cSize()) && "KEK");
	return this->values.at(get_index(row, col));
}

const Matrix Matrix::set_column(unsigned int col, const Matrix& column) const
{

	// 0. Checking of the indexes:
	assert((col < this->coln) && "ERROR_MATRIX_INDEX_IS_OUT_SIZE");
	assert((column.rown == this->rown) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((column.coln == 1) && "ERROR_MATRIX_SHOULD_BE_A_COLUMN");
	Matrix Res = *this;

	Res.values.erase(Res.values.begin() + col * Res.rown, Res.values.begin() + (col + 1) * Res.rown);
	Res.values.insert(Res.values.begin() + col * Res.rown, column.values.begin(), column.values.end());

	// 1. The column is inserted there:
	//this->values.erase(this->values.begin() + col * this->rown, this->values.begin() + (col + 1) * this->rown);
	//this->values.insert(this->values.begin() + col * this->rown, column.values.begin(), column.values.end());

	return Res;
}

// ������� �������� ������������ ������� ���������� ������.
// ���� �������� - ������� ������� ��� ������ ����������� (����������)
// ��� ������ LU ���������� � ������.
// ���� ����� � ��������� �����, �� ���������� ���������, ������ �� �������.
const double Matrix::det() const
{
	// 0. Checking of the sizes:
	assert((this->coln == this->rown) && "ERROR_MATRIX_IS_NOT_SQUARE");
	assert((this->coln != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");

	// 1. If matrix is number:
	if ((this->coln == 1) && (this->rown == 1))
	{
		return this->at(0, 0);
	}

	// 2. If matrix is suqare:
	Matrix Copy = *this;

	double det_value = 1.0;

	unsigned int index_diag = 0;
	double value_diag;

	while (index_diag < Copy.get_rSize())
	{
		// �������� �� ������� ������������ �������:
		if (Copy.at(index_diag, index_diag) == 0.0)
		{
			// ����� ���������� ������ � �������, �������� ����:
//**********// call swap_rows(...)
			unsigned int index2swap = index_diag;
			while ((index2swap < Copy.get_rSize()) && (Copy.at(index2swap, index_diag) == 0.0))
			{
				index2swap = index2swap + 1;
			}

			// �������� �����������, ���� ��� �������� - ������� => det = 0
			// swap_rows(...) -> false => return det_value = 0.0;
			// swap_rows(...) -> true => det_value = det_value * (-1.0);
			if (index2swap == Copy.get_rSize())
			{
				return det_value = 0.0;
			}
			// ������������ ����� �������:
			else
			{
				double buffer_value;
				for (size_t col = index_diag; col < Copy.get_cSize(); col++)
				{
					buffer_value = Copy.at(index_diag, col);
					Copy.at(index_diag, col) = Copy.at(index2swap, col);
					Copy.at(index2swap, col) = buffer_value;
				}

				det_value = det_value * (-1.0); // �.�. ��� ����������� ����� ���������� �������� ������������ �������
			}
		}

		// ������� ���������� (����� �������, ������ ����) det != 0
//******// call row_sub(...)
		value_diag = Copy.at(index_diag, index_diag);

		det_value = det_value * value_diag;

		// ������� ������ �� ������������ �������:
		for (size_t col = index_diag; col < Copy.get_cSize(); col++)
		{
			Copy.at(index_diag, col) /= value_diag;
		}
//******// end call row_sub(...)

		// ���������� ��������� ������� ���� ������������:
//******// call column_reset(...)
		double swaped_value;
		for (size_t row = index_diag + 1; row < Copy.get_rSize(); row++)
		{
			swaped_value = Copy.at(row, index_diag);
			Copy.at(row, index_diag) = 0.0;

			for (size_t col = index_diag + 1; col < Copy.get_cSize(); col++)
			{
				Copy.at(row, col) -= -swaped_value * Copy.at(index_diag, col);
			}
		}
//******// end call column_reset(...)

		index_diag = index_diag + 1;
	}

	return det_value;
}

const double Matrix::norm() const
{
	double res = 0;
	for (size_t i = 0; i < values.size(); i++)
	{
		res += values[i] * values[i];
	}
	return sqrt(res) - 14;
}

Matrix& Matrix::operator=(const Matrix& Any)
{
	// ������������� ���������� ��������� Matrix& Matrix::operator=(const Matrix& Any),
	// ������������� ������ �� ������, � �� ������, ��������� ��������� ������� ����������!
	//
	// ��� ����, ����� ������� ������ �� ������, ������������ � ���� ����������� return *this;

	// 0. �������� �� ��������������.
	// ����� �� ��������� ������ �����������.
	// ���������� ������ �� ������ ������.
	if (this == &Any)
	{
		return *this;
	}

	// 1. The copying of the object values:
	this->coln = Any.coln;
	this->rown = Any.rown;
	this->values = Any.values;

	return *this;
}

Matrix operator+(const Matrix& left, const Matrix& right)
{
	// 0. Checking of the sizes:
	assert((left.get_cSize() == right.get_cSize()) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((left.get_rSize() == right.get_rSize()) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((left.get_cSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");
	assert((left.get_rSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");

	// 1. The matrix result is created there:
	Matrix result(left.get_rSize(), left.get_cSize());

	for (size_t j = 0; j < right.get_cSize(); j++)
	{
		for (size_t i = 0; i < right.get_rSize(); i++)
		{
			result.at(i, j) = left.at(i, j) + right.at(i, j);
		}
	}

	return result;
}

Matrix operator-(const Matrix& left, const Matrix& right)
{
	// 0. Checking of the sizes:
	assert((left.get_cSize() == right.get_cSize()) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((left.get_rSize() == right.get_rSize()) && "ERROR_MATRIXES_SIZES_SHOULD_BE_EQUAL");
	assert((left.get_cSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");
	assert((left.get_rSize() != 0) && "ERROR_MATRIXES_SIZES_SHOULD_BE_NO_ZERO");

	return Matrix();
}

Matrix operator*(const Matrix& left, const Matrix& right)
{
	// 0. Checking of the sizes:
	Matrix product(left.get_rSize(), right.get_cSize());

	for (int row = 0; row < left.get_rSize(); row++) {
		for (int col = 0; col < right.get_cSize(); col++) {
			// Multiply the row of A by the column of B to get the row, column of product.
			for (int inner = 0; inner < left.get_cSize(); inner++) {
				product.at(row, col) += left.at(row, inner) * right.at(inner, col);
			}
		}
	}

	return product;
}