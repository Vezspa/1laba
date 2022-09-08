#include "Solver.h"
#include <iomanip>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>

void print(const Matrix& Any, unsigned int precicion)
{
	if ((Any.get_rSize() == 0) || (Any.get_cSize() == 0))
	{
		std::cout << "WARNING: printed matrix is empty!" << std::endl;
	}

	for (size_t i = 0; i < Any.get_rSize(); i++)
	{
		for (size_t j = 0; j < Any.get_cSize(); j++)
		{
			std::cout << std::setprecision(precicion) << std::scientific << Any.at(i, j) << "		";
		}
		std::cout << std::endl;
	}
}

Solver::Solver(std::string str) {
	read(str);
	SetB();
	LU = Decomposition(A);
}

Matrix Solver::Solve() {
	std::cout << "Methods: " << std::endl;
	std::cout << "1. Creamver" << std::endl;
	std::cout << "2. Decompose " << std::endl;

	int method;
	std::cin >> method;

	switch (method)
	{
	case 1:
		return Cramver();
		break;
	case 2:
		return Decompose();
		break;
	default:
		return Matrix();
		break;
	}
}

void Solver::read(std::string fullway2data)
{
	std::ifstream inputfile;
	inputfile.open(fullway2data);

	Matrix Res;

	if (inputfile.is_open())
	{
		std::string buff_s;
		double buff_d;
		std::vector <std::vector<double>> buff_data;
		std::vector <double> buff_data_row;

		while (getline(inputfile, buff_s))
		{
			std::istringstream buff_ss(buff_s);

			while (buff_ss >> buff_d)
			{
				buff_data_row.push_back(buff_d);
			}

			buff_data.push_back(buff_data_row);
			buff_data_row.clear();
		}

		Res = Matrix(buff_data.size(), buff_data.at(0).size());

		for (size_t row = 0; row < Res.get_rSize(); row++)
		{
			assert((buff_data.at(row).size() == Res.get_cSize()) && "ERROR_COPIED_MATRIX_COLUMNS_SIZES_SHOULD_BE_EQUAL");

			if (buff_data.at(row).size() != Res.get_cSize())
			{
				std::cout << "ERROR: copying matrix is failed! Process was stopped!" << std::endl;

				A = Res;
			}

			for (size_t col = 0; col < Res.get_cSize(); col++)
			{
				Res.at(row, col) = buff_data.at(row).at(col);
			}
		}
	}
	else
	{
		std::cout << "ERROR: copying matrix is failed! File isn't opened!" << std::endl;
	}

	A = Res;
}

void Solver::SetB() {
	Matrix x_exist = get_x_exist();

	printf("x_exist =\n");
	print(x_exist, 6);

	b = A * x_exist;
}

Matrix Solver::get_x_exist() {
	int rown = A.get_cSize();
	Matrix x_exist(rown, 1);

	for (size_t row = 1; row < rown; row++)
	{
		x_exist.at(row, 0) = std::sin(row + 1);
	}

	return x_exist;
}

Matrix Solver::Cramver() {

	Matrix Newx(A.get_rSize(), 1);
	const double det = A.det();
	if (det != 0) {
		std::vector <double> values;
		for (size_t i = 0; i < A.get_cSize(); i++)
		{
			Matrix insult = A.set_column(i, b);
			values.push_back(insult.det());
		}
		for (size_t i = 0; i < Newx.get_rSize(); i++)
		{
			Newx.at(i, 0) = values[i] / det;
		}
	}

	std::cout << std::endl << "Creamver method:" << std::endl;
	print(Newx, 6);

	return Newx;
}

const Matrix Solver::Decompose() const {

	int size = LU.get_size();

	Matrix y(size, 1);
	y.at(0, 0) = b.at(0, 0);

	for (size_t k = 1; k < size; k++)
	{
		y.at(k, 0) = b.at(k, 0);
		for (size_t p = 0; p < k; p++)
		{
			y.at(k, 0) -= LU.get_elemL(k, p) * y.at(p, 0);
		}
	}

	Matrix x(size, 1);
	int m = size - 1;
	x.at(m, 0) = y.at(m, 0) / LU.get_elemU(m, m);

	for (size_t k = size - 2; k >= 0 && k < size; k--)
	{
		x.at(k, 0) = y.at(k, 0);
		for (size_t p = k + 1; p < m; p++)
		{
			x.at(k, 0) -= LU.get_elemU(k, p) * x.at(p, 0);
		}
		x.at(k, 0) *= 1 / LU.get_elemU(k, k);
	}

	std::cout << std::endl << "LU decomposition method:" << std::endl;
	print(x, 6);

	return x;
}
