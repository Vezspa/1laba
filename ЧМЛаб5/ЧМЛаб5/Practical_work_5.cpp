#include <iostream>
#include <vector>
#include "Matrix.h"
#include "Solver.h"
#include <string>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>

int main()
{
	Solver solver("C:\\Program Files (x86)\\Работы\\Численные методы\\3.txt");

	Matrix first = solver.Solve();

	Matrix second = solver.Solve();

	std::cout << std::endl << "Method difference: " << (first - second).norm() << std::endl;

    return 0;
}