// Copyright 2023 Pashina Alina
#ifndef MODULES_TASK_4_PASHINA_A_SPARSE_MATRIX_STD_CRSMATRIX_STD_H_
#define MODULES_TASK_4_PASHINA_A_SPARSE_MATRIX_STD_CRSMATRIX_STD_H_

#include <string>
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>
#include "../../../3rdparty/unapproved/unapproved.h"

class CRSMatrix {
 public:
  int numRow, numCol;
  std::vector<double> valueCRS;
  std::vector<int> colsCRS;
  std::vector<int> pointerCRS;

  CRSMatrix(int numC, int numR, const std::vector<double>& myVal,
            const std::vector<int>& myColu, const std::vector<int>& myPointer);
  explicit CRSMatrix(int numC = 0, int numR = 0);
  explicit CRSMatrix(std::vector<std::vector<double>> matr);
  bool operator==(const CRSMatrix& matr) const;
  CRSMatrix MatrixTransp();
  CRSMatrix MatrixMult_threading(CRSMatrix m2);
};

std::vector<std::vector<double>> fillZero(int cols, int rows);
std::vector<std::vector<double>> createRandomMatrix(int cols, int rows,
                                                    double perc);
std::vector<std::vector<double>> multiplyVecMatrix(
    std::vector<std::vector<double>> myFirstMatrix,
    std::vector<std::vector<double>> mySecondMatrix);
void multiply_threading(CRSMatrix& result, const int thread_number,
                        const CRSMatrix& m1, const CRSMatrix& m2);

#endif  // MODULES_TASK_4_PASHINA_A_SPARSE_MATRIX_STD_CRSMATRIX_STD_H_
