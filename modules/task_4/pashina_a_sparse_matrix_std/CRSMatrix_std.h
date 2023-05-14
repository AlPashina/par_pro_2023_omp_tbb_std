// Copyright 2023 Pashina Alina
#ifndef MODULES_TASK_4_PASHINA_A_SPARSE_MATRIX_STD_CRSMATRIX_STD_H_
#define MODULES_TASK_4_PASHINA_A_SPARSE_MATRIX_STD_CRSMATRIX_STD_H_

#include <string>
#include <vector>

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

#endif  // MODULES_TASK_4_PASHINA_A_SPARSE_MATRIX_STD_CRSMATRIX_STD_H_
