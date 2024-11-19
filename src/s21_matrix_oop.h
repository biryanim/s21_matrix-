#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  int rows_, cols_;
  double **matrix_;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);
  ~S21Matrix();

  void SetRows(int rows);
  void SetCols(int cols);
  int GetRows();
  int GetCols();

  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  S21Matrix &operator=(const S21Matrix &right);
  friend const S21Matrix operator+(const S21Matrix &left,
                                   const S21Matrix &right);
  friend const S21Matrix operator-(const S21Matrix &left,
                                   const S21Matrix &right);
  friend const S21Matrix operator*(const S21Matrix &left,
                                   const S21Matrix &right);
  friend const S21Matrix operator*(const S21Matrix &left, const double num);
  friend const S21Matrix operator*(const double num, const S21Matrix &right);

  bool operator==(const S21Matrix &other);
  S21Matrix operator+=(const S21Matrix &other);
  S21Matrix operator-=(const S21Matrix &other);
  S21Matrix operator*=(const S21Matrix &other);
  S21Matrix operator*=(const double num);
  double &operator()(int i, int j);

  void Allocate(double ***matrix, int rows, int cols);
  void DeleteMatrix();
  int find_nonzero_rows(double **matrix, int rows, int columns);
  void swap_rows(double **matrix, int a, int b);
  S21Matrix ExpungeRowColumn(int row, int column);
  void SetSigns();
};
#endif