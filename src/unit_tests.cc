#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(Basic, defaultConstructor) {
  S21Matrix m;
  EXPECT_EQ(m.GetRows(), 1);
  EXPECT_EQ(m.GetCols(), 1);
}

TEST(Basic, parameterizedCconstructor) {
  S21Matrix m(2, 3);
  EXPECT_EQ(m.GetRows(), 2);
  EXPECT_EQ(m.GetCols(), 3);
}

TEST(Basic, moveConstructor) {
  S21Matrix m(2, 3);
  S21Matrix res(std::move(m));
  EXPECT_EQ(res.GetRows(), 2);
  EXPECT_EQ(res.GetCols(), 3);
  EXPECT_EQ(m.GetRows(), 0);
  EXPECT_EQ(m.GetCols(), 0);
}

TEST(GetterAndSetter, SetRows) {
  S21Matrix m(2, 3);
  m(1, 1) = 4.4;
  EXPECT_EQ(m(1, 1), 4.4);
  EXPECT_EQ(m.GetRows(), 2);
  EXPECT_EQ(m.GetCols(), 3);
  m.SetRows(4);
  m(3, 2) = 5.5;
  EXPECT_EQ(m(1, 1), 4.4);
  EXPECT_EQ(m(3, 2), 5.5);
  EXPECT_EQ(m.GetRows(), 4);
  EXPECT_EQ(m.GetCols(), 3);
}

TEST(GetterAndSetter, SetCols) {
  S21Matrix m(2, 3);
  m(1, 1) = 4.4;
  EXPECT_EQ(m(1, 1), 4.4);
  EXPECT_EQ(m.GetRows(), 2);
  EXPECT_EQ(m.GetCols(), 3);
  m.SetCols(5);
  m(1, 4) = 5.5;
  EXPECT_EQ(m(1, 1), 4.4);
  EXPECT_EQ(m(1, 4), 5.5);
  EXPECT_EQ(m.GetRows(), 2);
  EXPECT_EQ(m.GetCols(), 5);
}

TEST(assignmentOperator, brakets) {
  S21Matrix m(2, 3);
  m(1, 1) = 3;
  EXPECT_EQ(m(1, 1), 3);
}

TEST(functionalTest, copy) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(1, 1) = 1.1;
  b(1, 1) = 2.2;
  a = b;
  EXPECT_DOUBLE_EQ(a(1, 1), 2.2);
}

TEST(functionalTest, Plus) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(1, 1) = 1.1;
  b(1, 1) = 2.2;
  S21Matrix res = a + b;
  EXPECT_DOUBLE_EQ(res(1, 1), 3.3);
}

TEST(functionalTest, Plus2) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(1, 1) = 1.1;
  b(1, 1) = 2.2;
  a += b;
  EXPECT_DOUBLE_EQ(a(1, 1), 3.3);
}

TEST(functionalTest, Plus3) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(1, 1) = 1.1;
  b(1, 1) = 2.2;
  a.SumMatrix(b);
  EXPECT_DOUBLE_EQ(a(1, 1), 3.3);
}

TEST(functionalTest, Minus) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(1, 1) = 1.1;
  b(1, 1) = 2.2;
  S21Matrix res = a - b;
  EXPECT_DOUBLE_EQ(res(1, 1), -1.1);
}

TEST(functionalTest, Minus2) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(1, 1) = 1.1;
  b(1, 1) = 2.2;
  a -= b;
  EXPECT_DOUBLE_EQ(a(1, 1), -1.1);
}

TEST(functionalTest, Minus3) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(1, 1) = 1.1;
  b(1, 1) = 2.2;
  a.SubMatrix(b);
  EXPECT_DOUBLE_EQ(a(1, 1), -1.1);
}

TEST(functionalTest, MultMatrix) {
  S21Matrix a(3, 2);
  S21Matrix b(2, 3);
  a(1, 1) = 1.1;
  b(1, 1) = 2;
  S21Matrix res = a * b;
  EXPECT_DOUBLE_EQ(res(1, 1), 2.2);
}

TEST(functionalTest, MultMatrix2) {
  S21Matrix a(3, 2);
  S21Matrix b(2, 3);
  a(1, 1) = 1.1;
  b(1, 1) = 2;
  a *= b;
  EXPECT_DOUBLE_EQ(a(1, 1), 2.2);
}

TEST(functionalTest, MultMatrix3) {
  S21Matrix a(3, 2);
  S21Matrix b(2, 3);
  a(1, 1) = 1.1;
  b(1, 1) = 2;
  a.MulMatrix(b);
  EXPECT_DOUBLE_EQ(a(1, 1), 2.2);
}

TEST(functionalTest, MultMatrixNum) {
  S21Matrix a(3, 2);
  a(1, 1) = 1.1;
  S21Matrix res = a * 2;
  EXPECT_DOUBLE_EQ(res(1, 1), 2.2);
}

TEST(functionalTest, MultMatrixNum2) {
  S21Matrix a(3, 2);
  a(1, 1) = 1.1;
  a *= 2;
  EXPECT_DOUBLE_EQ(a(1, 1), 2.2);
}

TEST(functionalTest, MultMatrixNum3) {
  S21Matrix a(3, 2);
  a(1, 1) = 1.1;
  a.MulNumber(2);
  EXPECT_DOUBLE_EQ(a(1, 1), 2.2);
}

TEST(functionalTest, MultMatrixNum4) {
  S21Matrix a(3, 2);
  a(1, 1) = 1.1;
  S21Matrix b = 2 * a;
  EXPECT_DOUBLE_EQ(b(1, 1), 2.2);
}

TEST(functionalTest, equal) {
  S21Matrix a(2, 3);
  S21Matrix b(2, 2);
  EXPECT_EQ(a == b, false);
  b.SetCols(3);
  a(1, 1) = 1.1;
  b(1, 1) = 1.1;
  EXPECT_EQ(a == b, true);
  b(1, 2) = 1.1;
  EXPECT_EQ(a == b, false);
}

TEST(functionalFuncTest, transpose) {
  S21Matrix m(2, 3);
  m(0, 0) = 1;
  m(0, 1) = 2;
  m(0, 2) = 3;
  m(1, 0) = 4;
  m(1, 1) = 5;
  m(1, 2) = 6;
  // 1 2 3
  // 4 5 6
  EXPECT_EQ(m(1, 1), 5);
  m = std::move(m.Transpose());
  // 1 4
  // 2 5
  // 3 6
  EXPECT_EQ(m(0, 0), 1);
  EXPECT_EQ(m(2, 1), 6);
}

TEST(functionalFuncTest, determinant) {
  S21Matrix m(4, 4);
  m(0, 0) = 9;
  m(0, 1) = 2;
  m(0, 2) = 2;
  m(0, 3) = 4;

  m(1, 0) = 3;
  m(1, 1) = 4;
  m(1, 2) = 4;
  m(1, 3) = 4;

  m(2, 0) = 4;
  m(2, 1) = 4;
  m(2, 2) = 9;
  m(2, 3) = 9;

  m(3, 0) = 1;
  m(3, 1) = 1;
  m(3, 2) = 5;
  m(3, 3) = 1;
  EXPECT_EQ(m.Determinant(), -578);
  S21Matrix m1(1, 1);
  m1(0, 0) = 10;
  EXPECT_EQ(m1.Determinant(), 10);
  S21Matrix m2(2, 2);
  m2(0, 0) = 1.1;
  m2(0, 1) = 3.5;
  m2(1, 0) = -2;
  m2(1, 1) = 4;
  EXPECT_DOUBLE_EQ(m2.Determinant(), 11.4);
}

TEST(functionalFuncTest, determinant2) {
  S21Matrix m(4, 4);
  m(0, 0) = 1;
  m(0, 1) = 2;
  m(0, 2) = 3;
  m(0, 3) = 4;

  m(1, 0) = 1;
  m(1, 1) = 2;
  m(1, 2) = 5;
  m(1, 3) = 7;

  m(2, 0) = 1;
  m(2, 1) = 0;
  m(2, 2) = 6;
  m(2, 3) = 8;

  m(3, 0) = 1;
  m(3, 1) = 0;
  m(3, 2) = 6;
  m(3, 3) = 6;
  EXPECT_EQ(m.Determinant(), -8);
}

TEST(functionalFuncTest, calcComplements) {
  S21Matrix m(3, 3);
  m(0, 0) = 1;
  m(0, 1) = 2;
  m(0, 2) = 3;

  m(1, 0) = 0;
  m(1, 1) = 4;
  m(1, 2) = 2;

  m(2, 0) = 5;
  m(2, 1) = 2;
  m(2, 2) = 1;

  m = m.CalcComplements();

  EXPECT_EQ(m(0, 0), 0);
  EXPECT_EQ(m(0, 1), 10);
  EXPECT_EQ(m(0, 2), -20);
  EXPECT_EQ(m(1, 0), 4);
  EXPECT_EQ(m(1, 1), -14);
  EXPECT_EQ(m(1, 2), 8);
  EXPECT_EQ(m(2, 0), -8);
  EXPECT_EQ(m(2, 1), -2);
  EXPECT_EQ(m(2, 2), 4);
}

TEST(functionalFuncTest, inverseMatrix) {
  S21Matrix m(3, 3);
  // var 1
  m(0, 0) = 4;
  m(0, 1) = -2;
  m(0, 2) = 1;
  m(1, 0) = 1;
  m(1, 1) = 6;
  m(1, 2) = -2;
  m(2, 0) = 1;
  m(2, 1) = 0;
  m(2, 2) = 0;

  m = m.InverseMatrix();

  EXPECT_EQ(m(0, 1), 0);
  EXPECT_EQ(m(0, 2), 1);
  EXPECT_EQ(m(1, 0), 1);
  EXPECT_EQ(m(2, 0), 3);
  EXPECT_EQ(m(2, 1), 1);
  EXPECT_EQ(m(2, 2), -13);
}

TEST(functionalFuncTest, inverseMatrixEx) {
  S21Matrix m(3, 3);

  //  determ = 0
  m(0, 0) = 1;
  m(0, 1) = 1;
  m(0, 2) = 3;
  m(1, 0) = 4;
  m(1, 1) = 4;
  m(1, 2) = 6;
  m(2, 0) = 4;
  m(2, 1) = 4;
  m(2, 2) = 9;
  EXPECT_EQ(m.Determinant(), 0);
}

TEST(functionalFuncTest, inverseMatrixEx2) {
  S21Matrix m(3, 3);
  m(0, 0) = 1;
  m(0, 1) = 4;
  m(0, 2) = 1;
  m(1, 0) = 3;
  m(1, 1) = 7;
  m(1, 2) = 2;
  m(2, 0) = 3;
  m(2, 1) = 2;
  m(2, 2) = 1;
  EXPECT_EQ(m.Determinant(), 0);
}

TEST(functionalFuncTest, calcComplements2) {
  S21Matrix m(1, 1);
  S21Matrix n = m.CalcComplements();
  EXPECT_EQ(n(0, 0), 1);
}

TEST(functionalFuncTest, braketEx2) {
  S21Matrix m(3, 3);
  m(1, 1) = 1;
  EXPECT_EQ(m(1, 1), 1);
}

TEST(ThrowTest, assert_1) { ASSERT_ANY_THROW(S21Matrix m(-1, -1)); }

TEST(ThrowTest, assert_2) {
  S21Matrix m(34, 400);
  ASSERT_ANY_THROW(m.SetRows(-2));
  ASSERT_ANY_THROW(m.SetCols(-3));
}

TEST(ThrowTest, assert_3) {
  S21Matrix m(2, 5), n(4, 3);
  ASSERT_ANY_THROW(m.SumMatrix(n));
  ASSERT_ANY_THROW(n.SubMatrix(m));
  ASSERT_ANY_THROW(n.MulMatrix(m));
}

TEST(ThrowTest, assert_4) {
  S21Matrix m(1, 1);
  ASSERT_ANY_THROW(m.InverseMatrix());
}

TEST(ThrowTest, assert_5) {
  S21Matrix m(1, 1);
  ASSERT_ANY_THROW(m(-1, 0));
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}