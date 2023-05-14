// Copyright 2023 Pashina Alina
#include "../../../modules/task_4/pashina_a_sparse_matrix_std/CRSMatrix_std.h"

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../../../3rdparty/unapproved/unapproved.h"

CRSMatrix::CRSMatrix(int numC, int numR, const std::vector<double>& myVal,
                     const std::vector<int>& myColu,
                     const std::vector<int>& myPointer)
    : valueCRS(myVal),
      numCol(numC),
      colsCRS(myColu),
      numRow(numR),
      pointerCRS(myPointer) {}

CRSMatrix::CRSMatrix(int numC, int numR) {
  numCol = numC;
  numRow = numR;
}

CRSMatrix::CRSMatrix(std::vector<std::vector<double>> matr) {
  int indexCounter = 0;
  numRow = matr.size();
  numCol = matr[0].size();
  pointerCRS.push_back(indexCounter);
  for (int r = 0; r < numRow; r++) {
    for (int c = 0; c < numCol; c++) {
      if (matr[r][c] != 0) {
        valueCRS.push_back(matr[r][c]);
        indexCounter++;
        colsCRS.push_back(c);
      }
    }
    pointerCRS.push_back(indexCounter);
  }
}

bool CRSMatrix::operator==(const CRSMatrix& matr) const {
  if ((valueCRS == matr.valueCRS) && (numCol == matr.numCol) &&
      (colsCRS == matr.colsCRS) && (numRow == matr.numRow) &&
      (pointerCRS == matr.pointerCRS)) {
    return true;
  }
  return false;
}

CRSMatrix CRSMatrix::MatrixTransp() {
  CRSMatrix matr;
  std::vector<std::vector<int>> locCVec(numCol);
  std::vector<std::vector<double>> locVecVal(numCol);
  matr.numCol = numRow;
  int elemCounter = 0;
  matr.numRow = numCol;

  for (int r = 0; r < numRow; r++) {
    for (int ind = pointerCRS[r]; ind < pointerCRS[r + 1]; ind++) {
      int colInd = colsCRS[ind];
      locCVec[colInd].push_back(r);
      locVecVal[colInd].push_back(valueCRS[ind]);
    }
  }
  matr.pointerCRS.push_back(elemCounter);
  for (int col = 0; col < numCol; col++) {
    for (size_t ktmp = 0; ktmp < locCVec[col].size(); ktmp++) {
      matr.colsCRS.push_back(locCVec[col][ktmp]);
      matr.valueCRS.push_back(locVecVal[col][ktmp]);
    }
    elemCounter += locCVec[col].size();
    matr.pointerCRS.push_back(elemCounter);
  }
  return matr;
}

void multiply_threading(CRSMatrix& result, const int thread_number,
                        const CRSMatrix& m1, const CRSMatrix& m2) {
  const int n_elements = (m1.numRow * m2.numCol);
  const int n_operations = n_elements / thread_number;
  const int rest_operations = n_elements % thread_number;
  int start, end;
  if (thread_number < rest_operations) {
    start = (n_operations + 1) * thread_number;
    end = start + n_operations + 1;
  } else {
    start = n_operations * thread_number + rest_operations;
    end = start + n_operations;
  }

  std::vector<int> finColumn, finPointer;
  std::vector<double> finValue;
  int nRowZero = 0;
  finPointer.push_back(nRowZero);
  int finNumRows = m1.numRow;
  int finNumCols = m2.numCol;
  for (int r1 = start; r1 < end; r1++) {
    nRowZero = 0;
    for (int r2 = 0; r2 < m2.numRow; r2++) {
      int firstStart = m1.pointerCRS[r1];
      int secondStart = m2.pointerCRS[r2];
      double localSum = 0;
      while (((m2.pointerCRS[r2 + 1] - 1) >= secondStart) &&
             ((m1.pointerCRS[r1 + 1] - 1) >= firstStart)) {
        if (m1.colsCRS[firstStart] == m2.colsCRS[secondStart]) {
          localSum += (m1.valueCRS[firstStart] * m2.valueCRS[secondStart]);
          firstStart = firstStart + 1;
          secondStart = secondStart + 1;
        } else {
          if (m1.colsCRS[firstStart] < m2.colsCRS[secondStart]) {
            firstStart = firstStart + 1;
          } else {
            secondStart = secondStart + 1;
          }
        }
      }
      if (localSum != 0) {
        nRowZero = nRowZero + 1;
        finColumn.push_back(r2);
        finValue.push_back(localSum);
      }
    }
    finPointer.push_back(nRowZero + finPointer[r1]);
    for (int i = finPointer[r1]; i < finPointer[r1 + 1]; i++) {
      result.colsCRS.push_back(finColumn[i]);
      result.valueCRS.push_back(finValue[i]);
    }
  }
}

CRSMatrix CRSMatrix::MatrixMult_threading(CRSMatrix m2) {
  CRSMatrix result;
  if (numCol != m2.numCol) {
    throw std::runtime_error("Wrong matrix sizes!Change rows \n");
  }
  auto thread_num = std::thread::hardware_concurrency();
  result.numRow = numRow;
  result.numCol = m2.numCol;
  std::vector<std::thread> threads(thread_num);
  for (int i = 0; i < thread_num; ++i) {
    threads[i] = std::thread(multiply_threading, std::ref(result), i,
                             std::ref(*this), std::ref(m2));
  }
  for (int i = 0; i < thread_num; ++i) {
    threads[i].join();
  }
  return result;
}

std::vector<std::vector<double>> fillZero(int cols, int rows) {
  std::vector<std::vector<double>> res(rows);
  for (int m = 0; m < rows; m++) {
    for (int n = 0; n < cols; n++) {
      res[m].push_back(0);
    }
  }
  return res;
}
std::vector<std::vector<double>> createRandomMatrix(int cols, int rows,
                                                    double perc) {
  if (perc < 0 || perc > 1) {
    throw std::runtime_error("Wrong density \n");
  }
  std::random_device mydev;
  std::vector<std::vector<double>> res = fillZero(cols, rows);
  std::mt19937 gen(mydev());
  std::uniform_real_distribution<double> genP{0.0, 1.0};
  std::uniform_real_distribution<double> genVal{0.0, 25.0};
  for (int ro = 0; ro < rows; ro++) {
    for (int col = 0; col < cols; col++) {
      if (genP(gen) <= perc) {
        res[ro][col] = genVal(gen);
      }
    }
  }
  return res;
}

std::vector<std::vector<double>> multiplyVecMatrix(
    std::vector<std::vector<double>> myFirstMatrix,
    std::vector<std::vector<double>> mySecondMatrix) {
  int colsNum = mySecondMatrix[0].size();
  int rowsNumb = myFirstMatrix.size();
  std::vector<std::vector<double>> myResCRS = fillZero(colsNum, rowsNumb);
  for (int rr = 0; rr < rowsNumb; rr++) {
    for (int cc = 0; cc < colsNum; cc++) {
      myResCRS[rr][cc] = 0;
      for (size_t k = 0; k < myFirstMatrix[0].size(); k++) {
        myResCRS[rr][cc] += myFirstMatrix[rr][k] * mySecondMatrix[k][cc];
      }
    }
  }
  return myResCRS;
}
