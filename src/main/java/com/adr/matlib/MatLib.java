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

  public static double roundDouble(double d, int decimalPoint) {
    double p = 1;
    for(int i = 0; i < decimalPoint; i++) {
      p = roundDouble(p * 10);
    }

    return roundDouble(d * p) / p;
  }

  public static double roundDouble(double d) {
    long var2 = Double.doubleToRawLongBits(d);
    long var4 = (var2 & 9218868437227405312L) >> 52;
    long var6 = 1074L - var4;
    if((var6 & -64L) == 0L) {
      long var8 = var2 & 4503599627370495L | 4503599627370496L;
      if(var2 < 0L) {
        var8 = -var8;
      }

      return (var8 >> (int)var6) + 1L >> 1;
    } else {
      return (long)d;
    }
  }

  public static double[][][] gaussJordanElimination(double[][] matrixA, double[][] matrixB)
      throws NonConformableMatrixException{
    if(matrixB[0].length > 1) {
      throw new NonConformableMatrixException("Matrix B has width greater than 1");
    }

    int E = 1;
    double[][] matrixC = concatenateMatrix(matrixA, matrixB);

    for(int i = 0; i < matrixC.length; i++) {
      int p = computePivot(matrixC, i);

      if (roundDouble(matrixC[p][i]) == 0) {
        E = 0;
        return partitionMatrix(matrixC, matrixC.length - 1);
      }

      if (p > i) {
        matrixC = swapRow(matrixC, i, p);
      }

      matrixC = divideRow(matrixC, i, matrixC[i][i]);

      for (int j = 0; j < matrixC.length; j++) {
        if (j != i) {
          matrixC = subtractRow(matrixC, i, j, matrixC[j][i]);
        }
      }
    }

    return partitionMatrix(matrixC, matrixC.length - 1);
  }

  public static double[][][] determinants(double[][] matrixA, double[][] matrixB)
      throws NonConformableMatrixException {
    if (matrixB[0].length > 1) {
      throw new NonConformableMatrixException("Matrix B has width greater than 1");
    }

    int E = 1;
    double[][] matrixC = concatenateMatrix(matrixA, matrixB);

    for (int i = 0; i < matrixC.length; i++) {
      int p = computePivot(matrixC, i);

      if (roundDouble(matrixC[p][i]) == 0) {
        E = 0;
        return partitionMatrix(matrixC, matrixC.length - 1);
      }

      if (p > i) {
        matrixC = swapRow(matrixC, i, p);
      }

      matrixC = divideRow(matrixC, i, matrixC[i][i]);

      for (int j = 0; j < matrixC.length; j++) {
        if (j != i) {
          matrixC = subtractRow(matrixC, i, j, matrixC[j][i] / matrixC[j][j]);
        }
      }
    }

    return partitionMatrix(matrixC, matrixC.length - 1);
  }

  public static double[][][] invertMatrix(double[][] matrixA) {
    int E = 1;
    double[][] matrixC = concatenateMatrix(matrixA, generateIdentityMatrix(matrixA.length));

    for(int i = 0; i < matrixC.length; i++) {
      int p = computePivot(matrixC, i);

      if (roundDouble(matrixC[p][i]) == 0) {
        E = 0;
        return partitionMatrix(matrixC, matrixC.length / 2 + 1);
      }

      if(p > i) {
        matrixC = swapRow(matrixC, i, p);
      }

      matrixC = divideRow(matrixC, i, matrixC[i][i]);

      for(int j = 0; j < matrixC.length; j++) {
        if(j != i) {
          matrixC = subtractRow(matrixC, i, j, matrixC[j][i]);
        }
      }
    }

    return partitionMatrix(matrixC, matrixC.length / 2 + 1);
  }

  public static double[][][] partitionMatrix(double[][] matrix, int sizeOfFirstPart) {
    double[][] part1 = new double[matrix.length][sizeOfFirstPart + 1];
    double[][] part2 = new double[matrix.length][matrix[0].length - sizeOfFirstPart - 1];

    for(int i = 0; i < matrix.length; i++) {
      for(int j = 0; j < matrix[i].length; j++) {
        System.out.print(matrix[i][j] + "\t\t");
      }
      System.out.println("");
    }

    for(int i = 0; i < part1.length; i++) {
      for(int j = 0; j < part1[0].length; j++) {
        part1[i][j] = roundDouble(matrix[i][j], 10);
      }
    }

    for(int i = 0; i < part2.length; i++) {
      for(int j = 0; j < part2[0].length; j++) {
        part2[i][j] = roundDouble(matrix[i][j + sizeOfFirstPart + 1], 10);
      }
    }

    double[][][] parts = new double[2][][];
    parts[0] = part1;
    parts[1] = part2;

    return parts;
  }

  public static double[][] subtractRow(double[][] matrix, int row1, int row2, double times) {
    double[][] copy = matrix.clone();

    for(int i = 0; i < copy[0].length; i++) {
      copy[row2][i] -= times * copy[row1][i];
    }

    return copy;
  }

  public static double[][] divideRow(double[][] matrix, int row, double with) {
    double[][] copy = matrix.clone();

    for(int i = 0; i < copy[row].length; i++) {
      copy[row][i] = copy[row][i] / with * 1.0;
    }

    return copy;
  }

  public static double[][] swapRow(double[][] matrix, int replace, int with) {
    double[][] copy = matrix.clone();

    double[] temp = matrix[replace];
    copy[replace] = copy[with];
    copy[with] = temp;

    return copy;
  }

  public static int computePivot(double[][] matrix, int column) {
    int pivot = column;

    for(int i = column; i < matrix.length; i++) {
      if (matrix[i][column] > matrix[pivot][column]) {
        pivot = i;
      }
    }

    return pivot;
  }

  public static double[][] concatenateMatrix(double[][] matrixA, double[][] matrixB) {
    double[][] matrixC = new double[matrixA.length][matrixA[0].length + matrixB[0].length];

    for(int i = 0; i < matrixA.length; i++) {
      System.arraycopy(matrixA[i], 0, matrixC[i], 0, matrixA[i].length);
    }

    for (int i = 0; i < matrixB.length; i++) {
      System.arraycopy(matrixB[i], 0, matrixC[i], matrixA[i].length, matrixB[i].length);
    }

    return matrixC;
  }
}
