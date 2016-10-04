package com.adr.matlib;

import com.adr.matlib.exception.NonConformableMatrixException;

public final class MatLib {

  public static double[][] multiplyMatrix(double[][] matrixA, double[][] matrixB)
      throws NonConformableMatrixException {
    if (matrixA[0].length != matrixB.length) {
      throw new NonConformableMatrixException();
    }

    double[][] result = new double[matrixA.length][matrixB[0].length];

    for(int i = 0; i < result.length; i++) {
      result[i] = dotProduct(matrixA[i], matrixB);
    }

    return result;
  }

  public static double[] dotProduct(double[] a, double[][] b) {
    double[] newMatrix = new double[b[0].length];

   for(int i = 0; i < b[0].length; i++) {
     double sum = 0;

     for(int j = 0; j < a.length; j++) {
       sum += a[j] * b[j][i];
     }

     newMatrix[i] = sum;
   }

    return newMatrix;
  }

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
  public static double[][] multipleByScalar(double k, double[][] matrix) {
    double[][] cloneMatrix = matrix.clone();

    for(int i = 0; i < matrix.length; i++) {
      for(int j = 0; j < matrix[0].length; j++) {
        cloneMatrix[i][j] = k * cloneMatrix[i][j];
      }
    }

    return cloneMatrix;
  }

  public static double roundDouble(double d, int decimalPoint) {
    double p = 1;
    for(int i = 0; i < decimalPoint; i++) {
      p = roundDouble(p * 10);
    }

    return roundDouble(d * p) / p;
  }

  /**
   * Implementation of OpenJDK's round method for double precision numbers
   * @param d     Numbers being rounded
   * @return      Value of the argument rounded to the nearest long value
   */
  public static double roundDouble(double d) {
    // Double precision width
    int precisionWidth = 53;
    // 2^10 - 1
    int precisionBias = 1023;
    long signIfBitMask = 4503599627370495L;

    long rawLongBits = Double.doubleToRawLongBits(d);
    long biasedExp = (rawLongBits & Long.MAX_VALUE) >> precisionWidth - 1;
    long shift = (precisionWidth - 2 + precisionBias) - biasedExp;
    if((shift & -64L) == 0L) {
      long r = rawLongBits &  signIfBitMask | signIfBitMask + 1;
      if(rawLongBits < 0L) {
        r = -r;
      }

      return (r >> (int)shift) + 1L >> 1L;
    } else {
      return (long) d;
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

      if (matrixC[p][i] == 0) {
        matrixC[p][i] = 0;
        E = 0;
        return partitionMatrix(matrixC, matrixC[0].length - 1);
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

    return partitionMatrix(matrixC, matrixC[0].length - 1);
  }

  public static double[][][] gaussianElimination(double[][] matrixA, double[][] matrixB)
      throws NonConformableMatrixException {
    int E = 1;
    double[][] matrixC = concatenateMatrix(matrixA, matrixB);

    if (matrixC[0].length != matrixC.length + 1) {
      throw new NonConformableMatrixException("Matrix B has width greater than 1");
    }

    // For j = 0 to matrixC.length, do
    for (int j = 0; j < matrixC.length; j++) {
      // Compute pivot index j <= p <= matrixC.length
      int p = computePivot(matrixC, j);

      // If matrixC[p][j] == 0, set E = 0 and exit
      if (matrixC[p][j] == 0) {
        E = 0;
        return partitionMatrix(matrixC, matrixC[0].length - 1);
      }

      // If p > j, then swap rows p and j
      if (p > j) {
        matrixC = swapRow(matrixC, j, p);
      }

      // For each j > j, subtract matrixC[j][j] / matrixC[j][j] times row j from row j
      for (int i = 0; i < matrixC.length; i++) {
        if (i > j) {
          matrixC = subtractRow(matrixC, j, i, matrixC[i][j] / matrixC[j][j]);
        }
      }
    }

    // Partition matrix as C = [D, e] where D is n x n and e is n X 1
    return partitionMatrix(matrixC, matrixC[0].length - 1);
  }

  public static double[] backSubstitution(double[][][] partition) {
    double[][] D = partition[0];
    double[][] e = partition[1];
    double[]   x = new double[D.length];

    for(int j = D.length - 1; j >= 0; j--) {
      double sum = 0;

      for(int i = j + 1; i < D[j].length; i++) {
        sum += (D[j][i] * x[i]);
      }

      x[j] =  (e[j][0] - sum) / D[j][j];
    }

    return x;
  }

  public static double calculateDeterminant(double[][] matrix) throws NonConformableMatrixException {
    int r = 0;
    double det;
    double[][] tempMatrix = new double[matrix.length][matrix[0].length];

    for (int i = 0; i < matrix.length; i++) {
      System.arraycopy(matrix[i], 0, tempMatrix[i], 0, matrix[i].length);
    }

    if (tempMatrix.length != tempMatrix[0].length) {
      throw new NonConformableMatrixException("Matrix not n x n");
    }

    for(int j = 0; j < tempMatrix.length; j++) {
      int p = computePivot(tempMatrix, j);

      if(tempMatrix[p][j] == 0) {
        det = 0;
        return det;
      }

      if(p > j) {
        tempMatrix = swapRow(tempMatrix, j, p);
        r += 1;
      }

      for (int i = 0; i < tempMatrix.length; i++) {
        if (i > j) {
          subtractRow(tempMatrix, j, i, tempMatrix[i][j] / tempMatrix[j][j]);
        }
      }
    }

    double sum = 1;

    for(int i = 0; i < tempMatrix.length; i++) {
      sum *= tempMatrix[i][i];
    }

    return toPower(-1, r) * sum;
  }

  public static double toPower(double number, int power) {
    double newNumber = 1.0;
    if(power < 0) {
      for(int i = abs(power); i > 0; i--) {
        newNumber = newNumber / number;
      }
    } else {
      for(int i = power; i > 0; i--) {
        newNumber *= number;
      }
    }

    return newNumber;
  }

  public static double[][][] invertMatrix(double[][] matrixA) {
    int E = 1;
    double[][] matrixC = concatenateMatrix(matrixA, generateIdentityMatrix(matrixA.length));

    for(int i = 0; i < matrixC.length; i++) {
      int p = computePivot(matrixC, i);

      if (matrixC[p][i] == 0) {
        E = 0;
        return partitionMatrix(matrixC, matrixC[0].length / 2);
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

    return partitionMatrix(matrixC, matrixC[0].length / 2);
  }

  public static double[][][] partitionMatrix(double[][] matrix, int sizeOfFirstPart) {
    double[][] part1 = new double[matrix.length][sizeOfFirstPart];
    double[][] part2 = new double[matrix.length][matrix[0].length - sizeOfFirstPart];

    for (int i = 0; i < part1.length; i++) {
      System.arraycopy(matrix[i], 0, part1[i], 0, part1[i].length);
    }

    for (int i = 0; i < part2.length; i++) {
      System.arraycopy(matrix[i], 0 + sizeOfFirstPart, part2[i], 0, part2[i].length);
    }

    double[][][] parts = new double[2][][];
    parts[0] = part1;
    parts[1] = part2;

    return parts;
  }

  public static double[][] subtractRow(double[][] matrix, int value, int from, double times) {
    double[][] copy = matrix.clone();

    for(int i = 0; i < copy[0].length; i++) {
      copy[from][i] -= times * copy[value][i];
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
      if (abs(matrix[i][column]) > abs(matrix[pivot][column])) {
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

  public static double[][] transposeMatrix(double[][] matrix) {
    double[][] temp = new double[matrix[0].length][matrix.length];

   for(int i = 0; i < temp.length; i++) {
     for(int j = 0; j < temp[i].length; j++) {
       temp[i][j] = matrix[j][i];
     }
   }

    return temp;
  }
}
