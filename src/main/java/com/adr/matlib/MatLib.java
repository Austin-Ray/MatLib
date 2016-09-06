package com.adr.matlib;

import com.adr.matlib.exception.NonConformableMatrixException;

public final class MatLib {

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
  public static double[][] addMatrix(double[][] matrixA, double[][] matrixB) throws NonConformableMatrixException {
    if(matrixA.length != matrixB.length) {
      throw new NonConformableMatrixException();
    }

    double[][] matrixC = new double[matrixA.length][matrixA[0].length];

    for (int i = 0; i < matrixC.length; i++) {
      for (int j = 0; j < matrixC[0].length; j++) {
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
  public static double[][] subtractMatrix(double[][] matrixA, double[][] matrixB) throws NonConformableMatrixException {
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
}
