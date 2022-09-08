#pragma once
//Пользовательский тип(класс);
//Определены методы :
//	уметь извлекать элементы матрицы;
//	уметь определять размер матицы
//	операция вычисления определителя det();
//	метод замены столбца;
//	вычислять норму;
//	вычислять разность;
//	вычислять произведение;
//	копировать матрицы.

#include <iostream>
#include <vector>

class Matrix
{
private:
	// -1) The private geter gets a linear index:
	unsigned int get_index(unsigned int row, unsigned int col) const;	// Почему разумно сделать этот гетер приватным?
	// Такой подход позволяет скрыть реализацию хранения в виде линейного вектора.
	// В случае, если способ хранения по каким либо причинам изменить, то интерфейст доступа к элементам класса
	// меняться не будет!

	// 0) Values:
	std::vector <double> values;	// вектор значений (матрица будет хранится как последовательность столбцов)
	unsigned int rown;				// количество строк
	unsigned int coln;				// количество столбцов

public:
	// 1) Сonstructors:
	Matrix();
	Matrix(unsigned int rown, unsigned int coln); // инициализирует матрицу по умолчанию заданных размеров (заполненную нулями)

	// 2) Destructior:
	~Matrix();

	// 3) Geters and seters:
	const unsigned int get_rSize() const;
	const unsigned int get_cSize() const;
	double& at(unsigned int row, unsigned int col);
	const double& at(unsigned int row, unsigned int col) const;
	const Matrix set_column(unsigned int col, const Matrix& column) const; // должен изменять текущий объект

	// 4) Other methods:
	const double det() const;
	const double norm() const;

	// 5) Operathors:
	// добавить перегрузку операторов.
	// = Выполнено (через метод класса)
	// + Выполнено (перегрузка через функцию)
	// -		   (перегрузка через функцию)
	// *		   (перегрузка через функцию)
	// Использование перегрузки опреаторов через функции позволяет защитить данные от изменения.
	// Это достигается путём помещения функций перегрузки операторов за пределы поля класса.

	Matrix& operator=(const Matrix& Any); // Отличие от Matrix operator=(const Matrix& Any); смотри в теле
	// Оператор = лучше перегружать через метод класса. Это уменьшает код, позволяя еспользовать неявный указатель this
	// на объект стоящий слева от оператора.
	// Операторы +=, -= и т.д. тоже лучше перегружать через методы класса.

};

Matrix operator+(const Matrix& left, const Matrix& right);	// Исключение из тела класса, защищает данные от изменения.
Matrix operator-(const Matrix& left, const Matrix& right);	// Для операторов +,-,/, если есть возможность, лучше использовать
Matrix operator*(const Matrix& left, const Matrix& right);	// перегрузку через функции, как в данном примере.