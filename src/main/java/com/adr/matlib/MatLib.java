package com.adr.matlib;

import com.adr.matlib.exception.NonConformableMatrixException;

import java.util.Arrays;

public final class MatLib {

  private static void f(int s) {
    for(int i = 0; i < 512; i++) {
      double sum = 0.0;

      for(int k = 1; k <= s; k++) {
        sum += Math.sin(2 * Math.PI * (2 * k - 1) * (i * 1.0 / 512)) / (2 * k - 1);
      }

      System.out.println(sum);
    }
  }

  private static void g(int s) {
    for(int i = 0; i < 512; i++) {
      double sum = 0.0;

      for(int k = 1; k <= s; k++) {
        sum += Math.sin(2 * Math.PI * (2 * k) * (i * 1.0 / 512)) / (2 * k);
      }

      System.out.println(sum);
    }
  }

  public static double[] sinSummation(double a, double f1, double f2, double c, int samples) {
    double[] sum = new double[samples];

    for(int i = 0; i < samples; i++) {
      sum[i] = sin(a, f1, c, i * 1.0 / samples) + sin(a, f2, c, i * 1.0 / samples);
    }

    return sum;
  }

  public static double[] sinProduct(double a, double f1, double f2, double c, int samples) {
    double[] product = new double[samples];

    for(int i = 0; i < samples; i++) {
      product[i] = sin(a, f1, c, i * 1.0 / samples) * sin(a, f2, c, i * 1.0 / samples);
    }

    return product;
  }

  public static double sin(double a, double f, double c, double t) {
    return a * Math.sin(2 * Math.PI * f * (t - c));
  }

  public static double[] lowPassFilter(double[] data, int n) {
    double[] filter = new double[data.length];

    for(int i = 0; i < n; i++) {
      filter[i] = 1;
      filter[filter.length - i - 1] = 1;
    }

    return applyFilter(data, filter);
  }

  public static double[] highPassFilter(double[] data, int n) {
    double[] filter = new double[data.length];
    Arrays.fill(filter, 1);

    for(int i = 0; i < n; i++) {
      filter[i] = 0;
    }

    return applyFilter(data, filter);
  }

  public static double[] bandPassFilter(double[] data, int lowerBound, int upperBound) {
    return applyFilter(data, generateBandFilter(data.length, lowerBound, upperBound, 0, 1));
  }

  public static double[] notchFilter(double[] data, int lowerBound, int upperBound) {
    return applyFilter(data, generateBandFilter(data.length, lowerBound, upperBound, 1, 0));
  }

  private static double[] generateBandFilter(int length, int lowerBound, int upperBound, int fillValue, int filterValue) {
    double[] filter = new double[length];
    Arrays.fill(filter, fillValue);

    for(int i = lowerBound; i <= upperBound; i++) {
      filter[i] = filterValue;
    }

    for(int i = lowerBound; i < upperBound; i++) {
      filter[filter.length - 1 - i] = filterValue;
    }

    return filter;
  }

  private static double[] applyFilter(double[] data, double[] filter) {
    double[] duplicate = new double[data.length];
    System.arraycopy(data, 0, duplicate, 0, duplicate.length);

    for (int i = 0; i < duplicate.length; i++) {
      duplicate[i] *= filter[i];
    }

    return duplicate;
  }

  public static double[] PSD(Complex[] z) {
    double[] psd = new double[z.length];
    Complex[] fft = fastFourierTransform(z,1);

    for (int i = 0; i < fft.length; i++) {
      psd[i] = fft[i].times(fft[i]).abs();
    }

    return psd;
  }

  public static double[] fftConvolution(double[] uArray, int p) {
    int n = uArray.length;

    Complex[] u = new Complex[n];
    Complex[] h = new Complex[n];

    for (int i = 0; i < uArray.length; i++) {
      u[i] = new Complex(uArray[i], 0);
    }

    for (int i = 0; i < h.length; i++) {
      if(i < p) {
        h[i] = new Complex(1.0 / p, 0);
      } else {
        h[i] = new Complex(0, 0);
      }
    }

    Complex[] uFFT = fastFourierTransform(u, 1);
    Complex[] hFFT = fastFourierTransform(h, 1);

    Complex[] con = new Complex[n];

    for (int i = 0; i < con.length; i++) {
      con[i] = hFFT[i].times(uFFT[i]);
    }

    con = fastFourierTransform(con, -1);

    double[] real = new double[n];

    for (int i = 0; i < real.length; i++) {
      real[i] = con[i].re();
    }

    return real;
  }

  public static double[] normalizedCrossCorrelatiton(double[] y, double[] x) {
    double[] correlation = crossCorrelation(y, x);
    double[] autoXCorr = autoCorrelation(x);
    double[] autoYCorr = autoCorrelation(y);
    double[] norm = new double[correlation.length];

    for (int i = 0; i < norm.length; i++) {
      norm[i] = correlation[i] / Math.sqrt(autoXCorr[0] * autoYCorr[0]);
    }

    return norm;
  }

  public static double[] autoCorrelation(double[] x) {
    return crossCorrelation(x, x);
  }

  public static double[] crossCorrelation(double[] y, double[] x) {
    int n = x.length;
    double[] r = new double[n];
    double[] y2 = new double[n];
    double[] x2 = new double[n];

    System.arraycopy(y, 0, y2, 0, y.length);
    System.arraycopy(x, 0, x2, 0, x2.length);

    for (int i = 0; i < n; i++) {
      r[i] = 0.0;

      for(int k = i; k < n; k++) {
        r[i] += x2[k] * y2[k - i];
      }
    }

    return r;
  }

  /**
   * Performs a fast fourier transform
   * @param original     N x 1 Vector of complex numbers
   * @param d     Direction vector. 1 for FFT, -1 for inverse FFT
   * @return      FFT or inverse FFT of Z
   */
  public static Complex[] fastFourierTransform(Complex[] original, int d) {
    Complex[] z = new Complex[original.length];
    System.arraycopy(original, 0, z, 0, z.length);

    int n = z.length;
    double theta = (-2 * Math.PI * d) / n;
    int r = n / 2;
    int k, m;
    Complex t;

    for (int i = 1; i <= n-1; i *= 2) {
      Complex w = new Complex(Math.cos(i * theta), Math.sin(i * theta));

      for (k = 0; k <= n-1; k += 2*r) {
        Complex u = new Complex(1,0);
        for (m = 0; m <= r - 1; m++) {
          t = z[k+m].minus(z[k+m+r]);
          z[k+m] = z[k+m].plus(z[k+m+r]);
          z[k+m+r] = t.times(u);
          u = u.times(w);
        }
      }
      r = r / 2;
    }

    for (int i = 0; i <= n-1; i++) {
      r = i;
      k = 0;

      for(m = 1; m <= n - 1;) {
        k = 2*k + (r % 2);
        r /= 2;
        m *= 2;
      }

      if (k > i) {
        t = z[i];
        z[i] = z[k];
        z[k] = t;
      }
    }

    if (d < 0) {
      for(int i = 0; i <= n - 1; i++) {
        z[i] = z[i].scale(1.0/n);
      }
    }

    return z;
  }

  /**
   * Multiple two matrices
   * @param matrixA                           Matrix A
   * @param matrixB                           Matrix B
   * @return                                  Product of AB
   * @throws NonConformableMatrixException    Invalid matrix sizes
   */
  public static double[][] multiplyMatrix(double[][] matrixA, double[][] matrixB)
      throws NonConformableMatrixException {
    if (matrixA[0].length != matrixB.length) {
      throw new NonConformableMatrixException(String.format("Matrix A size: %s x %s, Matrix B size: %s x %s",
          matrixA.length, matrixA[0].length, matrixB.length, matrixB[0].length));
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
    double[][] cloneMatrix = copy2DMatrix(matrix);

    for(int i = 0; i < matrix.length; i++) {
      for(int j = 0; j < matrix[0].length; j++) {
        cloneMatrix[i][j] = k * cloneMatrix[i][j];
      }
    }

    return cloneMatrix;
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
    double[][] tempMatrix = copy2DMatrix(matrix);

    if (tempMatrix.length != tempMatrix[0].length) {
      throw new NonConformableMatrixException("Matrix not n x n");
    }

    if(matrix.length == 2) {
      return tempMatrix[0][0] * tempMatrix[1][1] - tempMatrix[0][1] * tempMatrix[1][0];
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
          tempMatrix = subtractRow(tempMatrix, j, i, tempMatrix[i][j] / tempMatrix[j][j]);
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
      for(int i = Math.abs(power); i > 0; i--) {
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
      System.arraycopy(matrix[i], sizeOfFirstPart, part2[i], 0, part2[i].length);
    }

    double[][][] parts = new double[2][][];
    parts[0] = part1;
    parts[1] = part2;

    return parts;
  }

  public static double[][] subtractRow(double[][] matrix, int value, int from, double times) {
    double[][] copy = copy2DMatrix(matrix);

    for(int i = 0; i < copy[0].length; i++) {
      copy[from][i] -= times * copy[value][i];
    }

    return copy;
  }

  public static double[][] divideRow(double[][] matrix, int row, double with) {
    double[][] copy = copy2DMatrix(matrix);

    for(int i = 0; i < copy[row].length; i++) {
      copy[row][i] = copy[row][i] / with * 1.0;
    }

    return copy;
  }

  public static double[][] swapRow(double[][] matrix, int replace, int with) {
    double[][] copy = copy2DMatrix(matrix);

    double[] temp = matrix[replace];
    copy[replace] = copy[with];
    copy[with] = temp;

    return copy;
  }

  public static int computePivot(double[][] matrix, int column) {
    int pivot = column;

    for(int i = column; i < matrix.length; i++) {
      if (Math.abs(matrix[i][column]) > Math.abs(matrix[pivot][column])) {
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

  /**
   * Transpose a matrix
   * @param matrix    Matrix being transposed
   * @return          Transpose of matrix
   */
  public static double[][] transposeMatrix(double[][] matrix) {
    double[][] temp = new double[matrix[0].length][matrix.length];

   for(int i = 0; i < temp.length; i++) {
     for(int j = 0; j < temp[i].length; j++) {
       temp[i][j] = matrix[j][i];
     }
   }

    return temp;
  }

  public static double traceMatrix(double[][] matrix) throws NonConformableMatrixException{
    if(matrix.length != matrix[0].length) {
      throw new NonConformableMatrixException("Not NxN");
    }

    double sum = 0;

    for(int i = 0; i < matrix.length; i++) {
      sum += matrix[i][i];
    }

    return sum;
  }

  public static double[][] eigenPowerMethod(double[][] matrix) throws NonConformableMatrixException {
    double[][] matrixA = copy2DMatrix(matrix);

    double e = 0.000001;
    double m = 100;
    double[][] y0 = generateVector(matrix.length, 1);
    double[][] r;
    double[][] eigenEst;
    double k = 0;

    double[][] x = multiplyMatrix(matrixA, y0);

    do {
      double[][] y = multipleByScalar(1.0 / vectorNorm(x), x);

      x = multiplyMatrix(matrixA, y);

      double[][] yTranpose = transposeMatrix(y);
      eigenEst = multiplyMatrix(multiplyMatrix(yTranpose, x), invertMatrix(multiplyMatrix(yTranpose, y))[1]);

      r = subtractMatrix(multipleByScalar(eigenEst[0][0], y), x);

      k+=1;
    }
    while (vectorNorm(r) > e && (k < m));

    return eigenEst;
  }

  public static double[] leverriersMethod(double[][] matrix) throws NonConformableMatrixException {
    int n = matrix.length - 1;

    // Matrix A
    double[][] tempMatrix = copy2DMatrix(matrix);

    // Create empty array for coefficients and B
    double[] coeffOfA = new double[tempMatrix.length];
    double[][][] bOfN = new double[tempMatrix.length][][];

    // Bn = A, where n is the length.
    bOfN[n] = copy2DMatrix(tempMatrix);

    // An = -trace(Bn)
    coeffOfA[n] = -1 * traceMatrix(bOfN[n]);

    for(int k = n - 1; k >= 0; k--) {
      // Get B(k+1)
      double[][] bOfK1 = copy2DMatrix(bOfN[k+1]);

      // Get a(k+1) and multiply the identity matrix
      double[][] aOfK1I = multipleByScalar(coeffOfA[k+1], generateIdentityMatrix(n + 1));

      // Add B(k+1) and the result of a(k+1) * I
      double[][] bOfK1PlusAOfK1I = addMatrix(bOfK1, aOfK1I);

      // A * (B(k+1) + a(k+1)*I)
      bOfN[k] = multiplyMatrix(tempMatrix, bOfK1PlusAOfK1I);

      double traceOfBk =  traceMatrix(bOfN[k]);
      double negateTraceOfBk = -1 * traceOfBk;
      double denominator = n - k + 1.0;

      coeffOfA[k] = negateTraceOfBk / denominator;
    }

    return coeffOfA;
  }

  public static double[][] findCovarianceMatrix(Vector[] classVectors, Vector mean, double scalar) throws NonConformableMatrixException {
    Vector[] tempVector = classVectors.clone();

    double[][] sumMatrix = new double[2][2];

    for(int i = 0; i < tempVector.length; i++) {
        double[][] v1 = MatLib.subtractMatrix(classVectors[i].matrix(), mean.matrix());
        double[][] transposeResult = MatLib.multiplyMatrix(v1, MatLib.transposeMatrix(v1));
        sumMatrix = MatLib.addMatrix(sumMatrix, transposeResult);
    }

    sumMatrix = MatLib.multipleByScalar(scalar, sumMatrix);

    return sumMatrix;
  }

  /**
   * Finds the L1 norm of a matrix
   * @param matrix    Matrix
   * @return          L1 norm of matrix
   */
  public static double matrixNorm(double[][] matrix) {
    double[][] tempMatrix = matrix.clone();
    double[] columnValues = new double[matrix[0].length];

    for (double[] aTempMatrix : tempMatrix) {
      for (int j = 0; j < aTempMatrix.length; j++) {
        columnValues[j] += Math.abs(aTempMatrix[j]);
      }
    }

    double max = Double.MIN_VALUE;

    for (double columnValue : columnValues) {
      if (columnValue > max) {
        max = columnValue;
      }
    }

    return max;
  }

  /**
   * Manual copy of a matrix
   * @param matrix    Matrix being copied
   * @return          Duplicate matrix
   */
  private static double[][] copy2DMatrix(double[][] matrix) {
    double[][] tempMatrix = new double[matrix.length][matrix[0].length];

    for (int i = 0; i < tempMatrix.length; i++) {
      System.arraycopy(matrix[i], 0, tempMatrix[i], 0, tempMatrix[i].length);
    }

    return tempMatrix;
  }

  public static double[][] normalizeVector(Vector u) {
    return multipleByScalar(1.0 / vectorNorm(u.matrix()), u.matrix());
  }

  /**
   * Euclidean norm for a vector
   * @param vector    Vector in matrix form
   * @return          Euclidean norm
   */
  public static double vectorNorm(double[][] vector) {
    double sum = 0.0;

    for (double[] aVector : vector) {
      sum += Math.pow(Math.abs(aVector[0]), 2);
    }

    return Math.sqrt(sum);
  }

  public static double[][][] jacobisMethod(double[][] matrixA) throws NonConformableMatrixException {
    if(matrixA.length != matrixA[0].length) {
      throw new NonConformableMatrixException("Not NxN");
    }

    int[] coords = highestMagnitudeCoordinate(matrixA);
    int p = coords[0];
    int q = coords[1];

    double[][] A = copy2DMatrix(matrixA);
    double[][] P = generateIdentityMatrix(matrixA.length);
    double[][] lambdas = new double[A.length][1];
    double e = 0.0001;

    while(Math.abs(A[p][q]) > e) {
      double phi = (1.0 / 2.0) * Math.atan((2 * A[p][q]) / (A[p][p] - A[q][q]));
      double[][] R = generateIdentityMatrix(A.length);

      // R(pp) = R(qq) = cos(phi)
      R[q][q] = Math.cos(phi);
      R[p][p] = R[q][q];

      // R(pq) = -R(qp) = -sin(phi)
      R[q][p] = -Math.sin(phi);
      R[p][q] = -R[q][p];

      P = multiplyMatrix(P, R);

      double[][] rTransponse = transposeMatrix(R);
      double[][] rTranposeA = multiplyMatrix(rTransponse, A);
      A = multiplyMatrix(rTranposeA, R);

    }

    for(int k = 0; k < lambdas.length; k++) {
      lambdas[k][0] = A[k][k];
    }

    double[][][] lambdaAndVectors = new double[2][][];
    lambdaAndVectors[0] = lambdas;
    lambdaAndVectors[1] = P;

    return lambdaAndVectors;
  }

  public static int[] highestMagnitudeCoordinate(double[][] matrix) {
    int[] coords = new int[2];
    double max = Double.MIN_VALUE;

    for(int row = 0; row < matrix.length; row++) {
      int maxColumnCoord = scanColumn(matrix, row);

      if(Math.abs(matrix[row][maxColumnCoord]) > max) {
        max = Math.abs(matrix[row][maxColumnCoord]);
        coords[0] = row;
        coords[1] = maxColumnCoord;
      }
    }

    return coords;
  }

  public static int scanColumn(double[][] matrix, int column) {
    double maxValue = Double.MIN_VALUE;
    int maxColumn = 0;

    for (int i = 0; i < matrix.length; i++) {
      if(Math.abs(matrix[i][column]) > maxValue) {
        maxValue = Math.abs(matrix[i][column]);
        maxColumn = i;
      }
    }

    return maxColumn;
  }

  public static double[][] generateVector(int length, int number) {
    double[][] vector = new double[length][1];

    for (int i = 0; i < vector.length; i++) {
      vector[i][0] = number;
    }

    return vector;
  }

  /**
   * Partitions a matrix to eigen vectors
   * @param matrix    Input matrix
   * @return          Array of Vectors
   */
  public static Vector[] matrixToEigenVectors(double[][] matrix) {
    Vector[] eigenVectors = new Vector[matrix.length];


    for(int column = 0; column < matrix[0].length; column++) {
      eigenVectors[column] = new Vector(matrix[0][column], matrix[1][column]);
    }

    return eigenVectors;
  }
}
