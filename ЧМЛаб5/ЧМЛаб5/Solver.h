#pragma once

#include "Matrix.h"
#include "Decomposition.h"

void print(const Matrix& Any, unsigned int precicion);

class Solver
{
private:
	Matrix A;
	Matrix b;

	Decomposition LU;

	Matrix get_x_exist();

	void read(std::string fullway2data);
	void SetB();

public:

	Solver(std::string str);

	Matrix Cramver();
	const Matrix Decompose() const;

	Matrix Solve();

};