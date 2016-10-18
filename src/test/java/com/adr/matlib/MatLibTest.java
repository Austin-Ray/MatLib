package com.adr.matlib;

import static org.junit.Assert.*;

import com.adr.matlib.exception.NonConformableMatrixException;
import org.junit.Test;

import java.util.Random;

public class MatLibTest {
  @Test
  public void highestMagnitudeCoordinate() throws Exception {
    double[][] matrixA = {{-1, 2}, {2, 2}};
    int[] expected = {0, 1};
    int[] result = MatLib.highestMagnitudeCoordinate(matrixA);

    checkArray(expected, result);
  }

  @Test
  public void scanColumn() throws Exception {
    double[][] matrix = {{1, 0}, {0, 1}};
    int expected1 = 0;

    assertEquals(expected1, MatLib.scanColumn(matrix, 0));
  }

  @Test
  public void jacobisMethod() throws Exception {
    double[][] matrixA = {{-1, 2},
                          {2, 2}};
    double[] expected = {3, -2};
    double[] result = MatLib.jacobisMethod(matrixA);

    //checkArray(expected, result, 0.001);
  }

  @Test
  public void normalizeVector() throws Exception {
    Vector vector = new Vector(8.0, -6.0);

    double[][] expected = {{4.0/5.0}, {-3.0/5.0}};
    double[][] result = MatLib.normalizeVector(vector);

    check2dArray(expected, result, 0.0001);
  }

  @Test
  public void traceMatrix() throws Exception {
    double[][] matrixA = {{1,1,1}, {1,1,1}, {1,1,1}};
    double expected = 3;

    assertEquals(expected, MatLib.traceMatrix(matrixA), 0.001);
  }

  @Test
  public void leverriersMethod() throws Exception {
    double[][] matrixA ={{2, -1, 1}, {-1, 2, 1}, {1, -1, 2}};
    double[] expected = {-6, 11, -6};
    double[] result = MatLib.leverriersMethod(matrixA);

    double[][] matrixB = {{1, -1, 0}, {0, 2, -1}, {-1, 0, 1}};
    double[] expected2 = {-1, 5, -4};
    double[] result2 = MatLib.leverriersMethod(matrixB);

    checkArray(expected, result, 0.0001);
    checkArray(expected2, result2, 0.001);
  }

  @Test
  public void matrixNorm() throws Exception {
    double[][] matrix = {{1, -7}, {-2, -3}};
    double[][] matrix2 = {{5, -4, 2}, {-1, 2, 3}, {-2, 1, 0}};
    double[][] matrix3 = {{1},{2}, {3}};

    double expected = 10;
    double expected2 = 8;
    double expected3 = 6;

    assertEquals(expected, MatLib.matrixNorm(matrix), 0.001);
    assertEquals(expected2, MatLib.matrixNorm(matrix2), 0.001);
    assertEquals(expected3, MatLib.matrixNorm(matrix3), 0.001);
  }

  @Test
  public void multiplyMatrix() throws Exception {
    double[][] matrixA = {{1, 2, 3},
                          {4, 5, 6}};
    double[][] matrixB = {{7, 8},
                          {9, 10},
                          {11, 12}};
    double[][] expected = {{58, 64}, {139, 154}};

    check2dArray(expected, MatLib.multiplyMatrix(matrixA, matrixB), 0.0001);
  }

  @Test
  public void transposeMatrix() throws Exception {
    double[][] matrix = {{1, 2}};
    double[][] expected = {{1}, {2}};

    assertArrayEquals(expected, MatLib.transposeMatrix(matrix));
  }

  @Test
  public void calculateDeterminant() throws Exception {
    double[][] matrix1 = {{3, 2},{5, 2}};
    assertEquals(-4, MatLib.calculateDeterminant(matrix1), 0.1);

    double[][] matrix = {{1, 4, 0}, {0, 2, 6}, {-1, 0, 1}};
    assertEquals(-22, MatLib.calculateDeterminant(matrix), 0.1);
  }

  @Test
  public void toPower() throws Exception {
    assertEquals(4, MatLib.toPower(2, 2), 0.1);
    assertEquals(1.0 / 4.0, MatLib.toPower(2, -2), 0.0001);
  }

  @Test
  public void backSubstitution() throws Exception {
    double[] expected = {-9, -2, 5};

    double[][] matrixA = {{1, 0, 2}, {2, -1, 3}, {4, 1, 8}};
    double[][] matrixB = {{1}, {-1}, {2}};

    double[][][] partitions = MatLib.gaussianElimination(matrixA, matrixB);

    double[] actual = MatLib.backSubstitution(partitions);

    assertArrayEquals(expected, actual, 0.001);
  }

  @Test
  public void gaussianElimination() throws Exception {
    double[][] matrixA = {{1, 0, 2}, {2, -1, 3}, {4, 1, 8}};
    double[][] matrixB = {{1}, {-1}, {2}};

    double[][][] partitions = MatLib.gaussianElimination(matrixA, matrixB);
    double[][] expectedArr = {{4, 1, 8}, { 0, -1.5, -1}, {0, 0, 1.0/6}};
    double[][] expectedResults = {{2}, {-2}, {5.0 / 6}};

    check2dArray(expectedArr, partitions[0], 0.01);
    check2dArray(expectedResults, partitions[1], 0.01);
  }

  @Test
  public void invertMatrix() throws Exception {
    double[][] matrixA = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
    double[][][] partitions = MatLib.invertMatrix(matrixA);

    double[][] expected = {{3.0/4, 1.0/2, 1.0/4}, {1.0/2, 1.0, 1.0/2}, {1.0/4, 1.0/2, 3.0/4}};
    double[][] expected2 = MatLib.generateIdentityMatrix(matrixA.length);

    check2dArray(expected2, partitions[0], 0.01);
    check2dArray(expected, partitions[1], 0.01);
  }

  @Test
  public void subtractRow() throws Exception {
    double[][] matrix = MatLib.generateIdentityMatrix(2);

    double[][] result = MatLib.subtractRow(matrix, 0, 1, 1);

    assertArrayEquals(new double[][]{{1, 0},{-1, 1}}, result);
  }

  @Test
  public void swapRow() throws Exception {
    double[][] matrix = MatLib.generateIdentityMatrix(2);
    matrix = MatLib.swapRow(matrix, 0, 1);
    assertArrayEquals(new double[][]{{0, 1}, {1, 0}}, matrix);
  }

  @Test
  public void computePivot() throws Exception {
    double[][] arr = MatLib.generateIdentityMatrix(2);
    MatLib.swapRow(arr, 0, 1);
    assertEquals(0, MatLib.computePivot(arr, 0));
  }

  @Test
  public void gaussJordanElimination() throws Exception {
    double[][] part1 = {{1, 1, 1},
                        {2, 3, 5},
                        {4, 0, 5}};

    double[][] part2 = {{5},
                        {8},
                        {2}};

    double[][][] partitions = MatLib.gaussJordanElimination(part1, part2);

    double[][] resultPart1 = partitions[0];
    double[][] resultPart2 = partitions[1];

    double[][] expected = { {1, 0, 0},
                            {0, 1, 0},
                            {0, 0, 1}};

    double[][] expected2 = new double[][]{{3}, {4}, {-2}};

    check2dArray(expected, resultPart1, 0.01);
    check2dArray(expected2, resultPart2, 0.01);
  }

  @Test
  public void concatenateMatrix() throws Exception {
    double[][] twoNGen = MatLib.generateIdentityMatrix(2);
    double[][] twoNGen2 = MatLib.generateIdentityMatrix(2);

    double[][] concatMatrix = MatLib.concatenateMatrix(twoNGen, twoNGen2);
    double[][] expectedMatrix = {{1,0,1,0},{0,1,0,1}};

    assertArrayEquals(expectedMatrix, concatMatrix);
  }

    @Test
  public void subtractMatrix() throws Exception {
    double[][] twoNGen = MatLib.generateIdentityMatrix(2);
    double[][] twoNGen2 =MatLib.generateIdentityMatrix(2);

    assertArrayEquals(new double[2][2], MatLib.subtractMatrix(twoNGen, twoNGen2));
  }

  @Test
  public void addMatrix() throws Exception {
    double[][] twoNGen = MatLib.generateIdentityMatrix(2);
    double[][] twoNAddExpected = {{2,0}, {0,2}};

    double[][] matrixA = {{5,0,2,3},{10,3,4,2}};
    double[][] matrixB = {{3,2,6,7},{10,2,5,7}};
    double[][] matrixC = {{8,2,8,10},{20,5,9,9}};

    assertArrayEquals(twoNAddExpected, MatLib.addMatrix(twoNGen, twoNGen));
    assertArrayEquals(matrixC, MatLib.addMatrix(matrixA, matrixB));

    double[][] mismatchMatrix = {{0}};
    double[][] mismatchMatrix2 = new double[2][6];

    try {
      MatLib.addMatrix(twoNGen, mismatchMatrix);
      MatLib.addMatrix(twoNGen, mismatchMatrix2);
      throw new Exception();
    } catch(NonConformableMatrixException e) {
      return;
    }
  }

  @Test
  public void multipleByScalar() throws Exception {
    double[][] twoNGen = MatLib.generateIdentityMatrix(2);
    double[][] fourNGen = MatLib.generateIdentityMatrix(4);

    double[][] twoNK = {{5,0},{0,5}};
    double[][] fourNK = {{8,0,0,0},{0,8,0,0},{0,0,8,0},{0,0,0,8}};

    assertArrayEquals(twoNK, MatLib.multipleByScalar(5, twoNGen));
    assertArrayEquals(fourNK, MatLib.multipleByScalar(8, fourNGen));
  }

  @Test
  public void generateIdentityMatrix() throws Exception {
    double[][] oneNGen = MatLib.generateIdentityMatrix(1);
    double[][] twoNGen = MatLib.generateIdentityMatrix(2);
    double[][] threeNGen = MatLib.generateIdentityMatrix(3);
    double[][] fourNGen = MatLib.generateIdentityMatrix(4);

    double[][] oneNConst = {{1}};
    double[][] twoNConst = {{1,0},{0,1}};
    double[][] threeNConst = {{1,0,0},{0,1,0},{0,0,1}};
    double[][] fourNConst = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

    assertArrayEquals(oneNConst, oneNGen);
    assertArrayEquals(twoNConst, twoNGen);
    assertArrayEquals(threeNConst, threeNGen);
    assertArrayEquals(fourNConst, fourNGen);
  }

  private void checkArray(double[] expected, double[] actual, double precision) {
    for (int i = 0; i < expected.length; i++) {
      assertEquals(expected[i], actual[i], precision);
    }
  }

  private void checkArray(int[] expected, int[] actual) {
    for (int i = 0; i < expected.length; i++) {
      assertEquals(expected[i], actual[i]);
    }
  }

  private void check2dArray(double[][] expected, double[][] actual, double precision) {
    for (int i = 0; i < expected.length; i++) {
      for (int j = 0; j < expected[i].length; j++) {
        assertEquals(expected[i][j], actual[i][j], precision);
      }
    }
  }
}