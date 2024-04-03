#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

class Matrix {
public:
int row1 = 0;
int cols1 = 0;
int row2 = 0;
int cols2 = 0;
std::vector<std::vector<double>> matrix1;
std::vector<std::vector<double>> matrix2; //вектора матриц
int N; //количество матриц

// Функция для вывода матрицы double
void printMatrix(const std::vector<std::vector<double>>& matrix) {
	for (const auto& row : matrix) {
		for (double element : row) {
			std::cout << element << " ";
		}
		std::cout << std::endl;
	}
}

// Функция для сложения
std::vector<std::vector<double>> AddMatrix(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) 
{
	if (N == 1) {
		row2 = row1;
		cols2 = cols1;
	}
	if ((row1 != row2) || (cols1 != cols2)) {
		std::cout << "size is not the same" << std::endl;
		exit(1);
	}
	std::vector<std::vector<double>> result;
	for(int i = 0; i < row1; i++) {
		std::vector<double> addmatrix;
		for(int j = 0; j < cols1; j++) {
			addmatrix.push_back(matrix1[i][j] + matrix2[i][j]);
		}
		result.push_back(addmatrix);
	}
	return result;
}

// Функция для вычитания
std::vector<std::vector<double>> SubMatrix(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) 
{
	if (N == 1) {
		row2 = row1;
		cols2 = cols1;
	}
	if ((row1 != row2) || (cols1 != cols2)) {
		std::cout << "size is not the same" << std::endl;
		exit(1);
	}
	std::vector<std::vector<double>> result;
	for(int i = 0; i < row1; i++) {
		std::vector<double> addmatrix;
		for(int j = 0; j < cols1; j++) {
			addmatrix.push_back(matrix1[i][j] - matrix2[i][j]);
		}
		result.push_back(addmatrix);
	}
	return result;
}

// Функция для умножения 
std::vector<std::vector<double>> MulMatrix(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) 
{
	if (N == 1) {
		row2 = row1;
		cols2 = cols1;
	}
	if (row1 != cols2) {
		std::cout << "the number of columns is not equal to the number of rows" << std::endl;
		exit(1);
	}
	std::vector<std::vector<double>> result;
	for(int i = 0; i < row1; i++) {
		std::vector<double> rez;
		for(int j = 0; j < cols2; j++) {
			double value = 0;
			for(int k = 0; k < cols1; k++) {
				value += matrix1[i][k] * matrix2[k][j];
			}
			rez.push_back(value);
		}
		result.push_back(rez);
	}
	return result;
} 

// Функция для деления
std::vector<std::vector<double>> DivMatrix(std::vector<std::vector<double>>& matrix1, std::vector<std::vector<double>>& matrix2) 
{
	if (N == 1) {
		row2 = row1;
		cols2 = cols1;
	}
	if (row1 != cols2) {
		std::cout << "the number of columns is not equal to the number of rows" << std::endl;
		exit(1);
	} 
	std::vector<std::vector<double>> result;
	std::vector<std::vector<double>> inversematrix = InverseMatrix(matrix2);
	for(int i = 0; i < row1; i++) {
		std::vector<double> rez;
		for(int j = 0; j < cols2; j++) {
			double value = 0;
			for(int k = 0; k < cols1; k++) {
				value += matrix1[i][k] * inversematrix[k][j];
			}
			rez.push_back(value);
		}
		result.push_back(rez);
	}
	return result;
} 

// Функция нахождения обратной матрицы
std::vector<std::vector<double>> InverseMatrix(std::vector<std::vector<double>>& matrix)
{
	if ((Determinant(matrix) == 0) or (row1 != cols1)) {
		std::cout << "determinant is null or the matrix is not square" << std::endl;
		exit(1);
	} else {
		std::vector<std::vector<double>> result;
		std::vector<std::vector<double>> matrix2 = TranspMatrix(matrix);
		for(int i = 0; i < row1; i++) {
			for(int j = 0; j < cols1; j++) {
				if (((i + 1) + (j + 1)) % 2 == 0) {
					matrix2[i][j] = matrix2[i][j];
				} else {
					matrix2[i][j] = -1 * matrix2[i][j];
				}
			}
		}
		std::vector<std::vector<double>> matrix3 = matrix2;
		for(int i = 0; i < row1; i++) {
			for(int j = 0; j < cols1; j++) {
				std::vector<std::vector<double>> res = GetMinor(matrix2, i, j);
				matrix3[i][j] = Determinant(res);
			}
		}
		matrix2 = matrix3;
		double count;
		count = 1 /  (double)Determinant(matrix);
		for(int i = 0; i < row1; i++) {
			std::vector<double> addmatrix;
			for(int j = 0; j < cols1; j++) {
				double c = matrix2[i][j] * count;
				addmatrix.push_back(c);
				
			}
			result.push_back(addmatrix);
		}
		return result;
	}
}


// Определитель (минор)
std::vector<std::vector<double>> GetMinor(std::vector<std::vector<double>>& matrix1, double row, double col) {
	std::vector<std::vector<double>> minor;
	int n = matrix1.size();
	for (int i = 0; i < n; ++i) {
		if (i == row) {
			continue;
		}
		std::vector<double> temp_row;
		for (int j = 0; j < n; ++j) {
			if (j != col) {
				temp_row.push_back(matrix1[i][j]);
			}
		}
		minor.push_back(temp_row);
	}
	return minor;
}
// Определитель
double Determinant(std::vector<std::vector<double>>& matrix1) {
	int n = matrix1.size();

	if (n == 1) {
		return matrix1[0][0];
	}

	double det = 0.0;
	double parity = 1.0;

	for (int i = 0; i < n; ++i) {
		std::vector<std::vector<double>> minor = GetMinor(matrix1, 0, i);
		det += parity * matrix1[0][i] * Determinant(minor);
		parity = -parity;
	}
	
	return det;
}

// Ранг матрицы
int RankMatrix(std::vector<std::vector<double>>& matrix) {
	int rank = 0;

	for (int i = 0; i < matrix.size(); ++i) {
		if (Determinant(matrix) != 0) {
			rank++;
		}
	}

	return rank;
}



// Транспонирование int
std::vector<std::vector<double>> TranspMatrix(std::vector<std::vector<double>>& matrix1) 
{
	std::vector<std::vector<double>> result(cols1,std::vector<double>(row1));
	for(int i = 0; i < row1; i++) {
		for(int j = 0; j < cols1; j++) {
			result[j][i] = matrix1[i][j];
		}
	}
	return result;
}

// Возведение матрицы в степень
std::vector<std::vector<double>> DegreeMatrix(std::vector<std::vector<double>>& matrix1, int degree) 
{
	if (row1 != cols1) {
		std::cout << "the number of columns is not equal to the number of rows" << std::endl;
		exit(1);
	} else {
		std::vector<std::vector<double>> result = matrix1;
		for(int i = 1; i < degree; i++) {
			result = MulMatrix(result, matrix1);
		}
		return result;
	}
}

// Определитель n-ого порядка
double DeterminantN(std::vector<std::vector<double>>& matrix1) {
	if (row1 != cols1) {
		std::cout << "the number of columns is not equal to the number of rows" << std::endl;
		exit(1);
	} else {
		int n = matrix1.size();
	
		if (n == 1) {
			return matrix1[0][0];
		}
	
		double det = 0.0;
		double parity = 1.0;
	
		for (int i = 0; i < n; ++i) {
			std::vector<std::vector<double>> minor = GetMinor(matrix1, 0, i);
			det += parity * matrix1[0][i] * Determinant(minor);
			parity = -parity;
		}
	
		return det;
	}
}


//функция считывания матриц
bool read_matrix(const std::string& filename) 
{
	std::ifstream file(filename);
	if (file.is_open()) {
		file >> N;
		if (N == 1) {
			file >> row1 >> cols1;
			for(int i = 0; i < row1; i++) {
				std::vector<double> matr;
				for(int j = 0; j < cols1; j++) {
					int value = 0;
					file >> value;
					matr.push_back(value);
				}
				matrix1.push_back(matr);
			}
		} else if (N == 2) {
			file >> row1 >> cols1 >> row2 >> cols2;
			for(int i = 0; i < row1; i++) {
				std::vector<double> matr;
				for(int j = 0; j < cols1; j++) {
					int value = 0;
					file >> value;
					matr.push_back(value);
				}
				matrix1.push_back(matr);
			}
			for(int i = 0; i < row2; i++) {
				std::vector<double> matr;
				for(int j = 0; j < cols2; j++) {
					int value = 0;
					file >> value;
					matr.push_back(value);
				}
				matrix2.push_back(matr);
			}
		} else {
			std::cout << "matrix is not 1 or 2" << std::endl;
			exit(1);
		}
		file.close();
		return true;
	} else {
		return false;
	}
}

};

int main() {
  	std::string  filename = "matrix.txt";
	Matrix matr;

	if (matr.read_matrix(filename)) {
		std::cout << "success" << std::endl;
	} else {
		std::cout << "not success" << std::endl;
	}
	//matr.printMatrix(matr.InverseMatrix(matr.matrix1));
	//std::cout << matr.Determinant(matr.matrix1) << std::endl;
	//std::cout << matr.RankMatrix(matr.matrix1) << std::endl;
	//matr.printMatrix(matr.DegreeMatrix(matr.matrix1, 2));
	//matr.printMatrix(matr.AddMatrix(matr.matrix1, matr.matrix1));
	//matr.printMatrix(matr.SubMatrix(matr.matrix1, matr.matrix1));
	matr.printMatrix(matr.InverseMatrix(matr.matrix1));
	matr.printMatrix(matr.DivMatrix(matr.matrix1, matr.matrix1));
	matr.printMatrix(matr.MulMatrix(matr.matrix1, matr.matrix1));
	return 0;
}