#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(1), cols_(1), matrix_(nullptr) {
  Allocate(&matrix_, rows_, cols_);
}

S21Matrix::S21Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), matrix_(nullptr) {
  Allocate(&matrix_, rows_, cols_);
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] = other.matrix_[i][j];
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() { DeleteMatrix(); }

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) return false;
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) return false;
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument("Different dimensionality of matrices.");
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] += other.matrix_[i][j];
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument("Different dimensionality of matrices.");
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] -= other.matrix_[i][j];
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] *= num;
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_)
    throw std::invalid_argument(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix.");
  S21Matrix res(rows_, other.cols_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < other.cols_; j++)
      for (int k = 0; k < cols_; k++)
        res.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
  *this = res;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix transpose(cols_, rows_);
  for (int i = 0; i < cols_; i++)
    for (int j = 0; j < rows_; j++) transpose.matrix_[i][j] = matrix_[j][i];
  return transpose;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) throw std::invalid_argument("The matrix is not square.");
  S21Matrix copy(*this);
  double res = 1;
  int sign = 1;
  for (int i = 0; i < rows_ && true; i++) {
    if (copy.matrix_[i][i] == 0) {
      int swap;
      if (!(swap = find_nonzero_rows(copy.matrix_, i, i))) {
        return 0;
      } else {
        swap_rows(copy.matrix_, i, swap);
        sign *= -1;
      }
    }
    for (int j = i + 1; j < rows_; j++) {
      double div = copy.matrix_[j][i] / copy.matrix_[i][i];
      if (div == 0) continue;
      for (int k = i; k < rows_; k++) {
        copy.matrix_[j][k] -= div * copy.matrix_[i][k];
      }
    }
  }
  for (int i = 0; i < rows_; i++) res *= copy.matrix_[i][i];
  res *= sign;
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) throw std::invalid_argument("The matrix is not square.");
  S21Matrix res(rows_, cols_);
  if (rows_ == 1)
    res.matrix_[0][0] = 1.;
  else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < rows_; j++) {
        res.matrix_[i][j] = ExpungeRowColumn(i, j).Determinant();
      }
    }
    res.SetSigns();
  }
  return res;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (fabs(det) < 1e-7)
    throw std::runtime_error("The determinant of the matrix is 0.");
  S21Matrix res(rows_, cols_);
  res = CalcComplements().Transpose();
  res.MulNumber(1 / det);
  return res;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &right) {
  DeleteMatrix();
  rows_ = right.rows_;
  cols_ = right.cols_;
  Allocate(&matrix_, rows_, cols_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) matrix_[i][j] = right.matrix_[i][j];
  return *this;
}

const S21Matrix operator+(const S21Matrix &left, const S21Matrix &right) {
  S21Matrix res(left);
  res.SumMatrix(right);
  return res;
}

const S21Matrix operator-(const S21Matrix &left, const S21Matrix &right) {
  S21Matrix res(left);
  res.SubMatrix(right);
  return res;
}

const S21Matrix operator*(const S21Matrix &left, const S21Matrix &right) {
  S21Matrix res(left);
  res.MulMatrix(right);
  return res;
}

const S21Matrix operator*(const S21Matrix &left, const double num) {
  S21Matrix res(left);
  res.MulNumber(num);
  return res;
}
const S21Matrix operator*(const double num, const S21Matrix &right) {
  S21Matrix res(right);
  res.MulNumber(num);
  return res;
}

bool S21Matrix::operator==(const S21Matrix &other) { return EqMatrix(other); }

S21Matrix S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || i < 0 || j >= cols_ || j < 0)
    throw std::invalid_argument("Index outside the matrix.");
  return matrix_[i][j];
}

void S21Matrix::SetCols(int cols) {
  if (cols_ == cols) return;
  double **tmp;
  Allocate(&tmp, rows_, cols);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_ && j < cols; j++) tmp[i][j] = matrix_[i][j];
  }
  DeleteMatrix();
  matrix_ = tmp;
  cols_ = cols;
}

void S21Matrix::SetRows(int rows) {
  if (rows_ == rows) return;
  double **tmp;
  Allocate(&tmp, rows, cols_);
  for (int i = 0; i < rows && i < rows_; i++) {
    for (int j = 0; j < cols_; j++) tmp[i][j] = matrix_[i][j];
  }
  DeleteMatrix();
  matrix_ = tmp;
  rows_ = rows;
}

int S21Matrix::GetRows() { return rows_; }

int S21Matrix::GetCols() { return cols_; }

void S21Matrix::SetSigns() {
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      matrix_[i][j] *= ((i + j) % 2 == 0) ? 1 : -1;
}

S21Matrix S21Matrix::ExpungeRowColumn(int row, int column) {
  S21Matrix tmp(rows_ - 1, cols_ - 1);
  for (int i = 0, m = 0; i < rows_; i++) {
    for (int j = 0, n = 0; j < cols_; j++) {
      if (i != row && j != column) {
        tmp.matrix_[m][n] = matrix_[i][j];
        n++;
        if (n >= tmp.cols_) {
          n = 0;
          m++;
        }
      }
    }
  }
  return tmp;
}

void S21Matrix::swap_rows(double **matrix, int a, int b) {
  double *tmp = matrix[a];
  matrix[a] = matrix[b];
  matrix[b] = tmp;
}

int S21Matrix::find_nonzero_rows(double **matrix, int rows, int columns) {
  bool flag = false;
  int ind = 0;
  for (int i = rows; i < rows_ && !flag; i++)
    if (matrix[i][columns] != 0) {
      flag = true;
      ind = i;
    }
  return ind;
}

void S21Matrix::Allocate(double ***matrix, int rows, int cols) {
  if (rows < 1 || cols < 1)
    throw std::invalid_argument("Invalid size of rows or columns");
  *matrix = new double *[rows]();
  for (int i = 0; i < rows; i++) (*matrix)[i] = new double[cols]();
}

void S21Matrix::DeleteMatrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) delete[] matrix_[i];
    delete[] matrix_;
  }
}
