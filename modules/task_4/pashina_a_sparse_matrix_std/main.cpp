// Copyright 2023 Pashina Alina
#include <gtest/gtest.h>

#include <vector>

#include "../../../3rdparty/unapproved/unapproved.h"
#include "./CRSMatrix_std.h"

TEST(MatrixCRS_support, make_and_fill_zero_matrix) {
  EXPECT_NO_THROW(fillZero(5, 5));
}
TEST(MatrixCRS_support, random_vector_matrix) {
  EXPECT_NO_THROW(createRandomMatrix(5, 5, 0.5));
}
TEST(MatrixCRS_support, vector_matrix_big_density) {
  EXPECT_ANY_THROW(createRandomMatrix(5, 5, 2));
}
TEST(MatrixCRS_support, vector_matrix_oppos_density) {
  EXPECT_ANY_THROW(createRandomMatrix(5, 5, -2));
}
TEST(CRSMatrix, matrix_mult_vec) {
  std::vector<std::vector<double>> ve1 = {{0, 8}, {0, 2}};
  std::vector<std::vector<double>> ve2 = {{1, 0}, {1, 0}};
  std::vector<std::vector<double>> ve3 = {{8, 0}, {2, 0}};
  EXPECT_EQ(multiplyVecMatrix(ve1, ve2), ve3);
}

TEST(CRSMatrix, create_crs_matrix) {
  std::vector<std::vector<double>> veMatrix = {
      {0, 5},
      {9, 0},
      {0, 6},
  };
  EXPECT_NO_THROW(CRSMatrix resCRS(veMatrix));
}

TEST(CRSMatrix, can_transp_CRSMatrix) {
  std::vector<double> row1 = {0, 0, 0, 8};
  std::vector<double> row2 = {6, 0, 0, 0};
  std::vector<double> row3 = {0, 15, 0, 0};
  std::vector<double> row4 = {6, 0, 0, 0};
  std::vector<double> row5 = {0, 1, 0, 72};
  std::vector<std::vector<double>> vT = {row1, row2, row3, row4, row5};

  std::vector<double> res1 = {0, 6, 0, 6, 0};
  std::vector<double> res2 = {0, 0, 15, 0, 1};
  std::vector<double> res3 = {0, 0, 0, 0, 0};
  std::vector<double> res4 = {8, 0, 0, 0, 72};
  std::vector<std::vector<double>> resMatrix = {res1, res2, res3, res4};
  EXPECT_EQ(CRSMatrix(vT).MatrixTransp(), CRSMatrix(resMatrix));
}
TEST(CRSMatrix, random_mult_CRSMatrix) {
  std::vector<std::vector<double>> veM1 = createRandomMatrix(3, 3, 0.4);
  std::vector<std::vector<double>> veM2 = createRandomMatrix(3, 3, 0.4);
  CRSMatrix myCRSMatr1(veM1);
  CRSMatrix myCRSMatr2(veM2);
  EXPECT_EQ(myCRSMatr1.MatrixMult_threading(myCRSMatr2.MatrixTransp()),
            CRSMatrix(multiplyVecMatrix(veM1, veM2)));
}
TEST(CRSMatrix, expected_mult_CRSMatrix) {
  CRSMatrix firstM(6, 6, {1, 2, 3, 4, 2, 5, 7}, {0, 4, 2, 3, 3, 5, 1},
                   {0, 2, 4, 4, 6, 6, 7});
  CRSMatrix secondM(6, 6, {3, 1, 1, 3, 1, 3, 2, 4}, {0, 5, 2, 3, 0, 2, 1, 4},
                    {0, 2, 3, 3, 4, 6, 8});
  CRSMatrix resM(6, 6, {5, 6, 1, 12, 10, 6, 20, 7}, {0, 2, 5, 3, 1, 3, 4, 2},
                 {0, 3, 4, 4, 7, 7, 8});
  CRSMatrix res = firstM.MatrixMult_threading(secondM.MatrixTransp());
  EXPECT_EQ(res, resM);
}
TEST(CRSMatrix, mult_CRSMatrix_bad_size) {
  std::vector<std::vector<double>> veM1 = createRandomMatrix(3, 2, 0.4);
  std::vector<std::vector<double>> veM2 = createRandomMatrix(5, 4, 0.4);
  CRSMatrix myCRSMatr1(veM1);
  CRSMatrix myCRSMatr2(veM2);
  EXPECT_ANY_THROW(myCRSMatr1.MatrixMult_threading(myCRSMatr2.MatrixTransp()));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
