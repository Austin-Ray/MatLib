package com.adr.matlib;

import com.adr.matlib.exception.NonConformableMatrixException;

final class MatLib {

  /**
   * Absolute value of an integer value
   * @param n     number
   * @return      Absolute value of a number
   */
  public static int abs(int n) {
    return n > 0 ? n : -n;
  }

  /**
   * Absolute value of a long number
   * @param n     number
   * @return      absolute value of number
   */
  public static long abs(long n) {
    return n > 0L ? n : n * -1;
  }

  /**
   * Absolute value of a floating point value with single precision
   * @param n     number
   * @return      Absolute value of number
   */
  public static float abs(float n) {
    return n >= 0.0F ? n : 0.0F - n;
  }

  /**
   * Absolute value of a floating point value with double precision
   * @param n     number
   * @return      absolute value of number
   */
  public static double abs(double n) {
    return n >= 0.0D ? n : 0.0D - n;
  }

  /**
   * Adds two matrices and returns the result
   * @param matrixA                             First matrix
   * @param matrixB                             Second matrix
   * @return                                     First matrix + second matrix
   * @throws NonConformableMatrixException       Matrix do not have compatible sizes
   */
  public static double[][] addMatrix(double[][] matrixA, double[][] matrixB)
      throws NonConformableMatrixException {
    if(matrixA.length != matrixB.length) {
      throw new NonConformableMatrixException();
    }

    double[][] matrixC = new double[matrixA.length][matrixA[0].length];

    for (int i = 0; i < matrixC.length; i++) {
      for (int j = 0; j < matrixC[i].length; j++) {
        if(matrixA[i].length != matrixB[i].length) {
          throw new NonConformableMatrixException();
        }

        matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
      }
    }

    return matrixC;
  }

  /**
   * Subtract Matrix B from Matrix A
   * @param matrixA                           Matrix A
   * @param matrixB                           Matrix B
   * @return                                  Resulting Matrix C
   * @throws NonConformableMatrixException    When Matrix sizes are not compatible.
   */
  public static double[][] subtractMatrix(double[][] matrixA, double[][] matrixB)
      throws NonConformableMatrixException {
    return addMatrix(matrixA, multipleByScalar(-1, matrixB));
  }

  /**
   * Generate a n x n identity matrix
   * @param n     Column/Row size
   * @return      n x n identity matrix
   */
  public static double[][] generateIdentityMatrix(int n) {
    double[][] identityMatrix = new double[n][n];

    for(int i = 0; i < n; i++) {
      identityMatrix[i][i] = 1;
    }

    return identityMatrix;
  }

  /**
   * Multiply a matrix by a scalar vector
   * @param k         Scalar multiple
   * @param matrix    Original matrix
   * @return          Scaled matrix
   */
  public static double[][] multipleByScalar(int k, double[][] matrix) {
    for(int i = 0; i < matrix.length; i++) {
      for(int j = 0; j < matrix[0].length; j++) {
        matrix[i][j] = k * matrix[i][j];
      }
    }

    return matrix;
  }

  /**
   * Gauss-Jordan elimination for a matrix system Ax=B
   * @param matrixA     Matrix A
   * @param matrixB     Matrix B
   * @return            Solution for a system of equations
   */
  public static double[][] gaussJordanElimination(double[][] matrixA, double[][] matrixB) {
    int E = 1;
    double[][] augMatrix = concatenateMatrix(matrixA, matrixB);

    for(int j = 0; j < augMatrix.length; j++) {
      int p = computePivot(augMatrix, j);

      if(augMatrix[p][j] == 0) {
        E = 0;
        return new double[][] {{-1}};
      }

      swapRow(augMatrix, j, p);
      divideRow(augMatrix, j);

      for(int i = 0; i < augMatrix.length; i++) {
        if(i != j) {
          subtractRow(augMatrix, j, i, augMatrix[j][i]);
        }
      }
    }

    return augMatrix;
  }

  public static void subtractRow(double[][] matrix, int row1, int row2, double times) {
    for(int i = 0; i < matrix[row1].length; i++) {
      matrix[row2][i] -= times * matrix[row1][i];
    }
  }

  /**
   * Divides a row in a matrix by a value
   * @param matrix    Matrix being operated on
   * @param row       Row being divided
   */
  public static void divideRow(double[][] matrix, int row) {
    for(int i = 0; i < matrix[row].length; i++ ) {
      matrix[row][i] /= matrix[row][row];
    }
  }

  /**
   * Swaps two rows in a matrix
   * @param matrix    Matrix operation is being performed on
   * @param row1      Row being replaced
   * @param row2      Row replacing the other
   */
  public static void swapRow(double[][] matrix, int row1, int row2) {
    double[] temp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = temp;
  }

  /**
   * Computes a pivot point by checks every row in a column for
   * the greatest value and returns the row number for that value
   * @param matrix    Matrix being checked
   * @param column    Column to check in
   * @return          Row # with the greatest value
   */
  public static int computePivot(double[][] matrix, int column) {
    int p = Integer.MIN_VALUE;

    for(int i = 0; i < matrix.length; i++) {
      if(matrix[i][column] > p) {
        p = i;
      }
    }

    return p;
  }

  public static double[][] concatenateMatrix(double[][] matrixA, double[][] matrixB) {
    double[][] matrixC = new double[matrixA.length][matrixA[0].length + matrixB[0].length];

    for(int i = 0; i < matrixA.length; i++) {
      for(int j = 0; j < matrixA[i].length; j++) {
        matrixC[i][j] = matrixA[i][j];
      }
    }

    for (int i = 0; i < matrixB.length; i++) {
      for(int j = 0; j < matrixB[i].length; j++) {
        matrixC[i][j + matrixA[i].length] = matrixB[i][j];
      }
    }

    return matrixC;
  }
}
