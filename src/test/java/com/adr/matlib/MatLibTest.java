package com.adr.matlib;

import static org.junit.Assert.*;

import com.adr.matlib.exception.NonConformableMatrixException;
import org.junit.Test;

import java.util.Random;

public class MatLibTest {
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
    double[][] mismatchMatrix2 = {{0, 0, 0}, {0}};

    try {
      MatLib.addMatrix(twoNGen, mismatchMatrix);
      MatLib.addMatrix(twoNGen, mismatchMatrix2);
      throw new Exception();
    } catch(NonConformableMatrixException e) {
      return;
    }
  }

  @Test
  public void abs() throws Exception {
    int i = Integer.MAX_VALUE;
    int j = Integer.MIN_VALUE;

    // Overflows because the min value has a larger
    // absolute value than the max value of a int
    assertEquals(Integer.MAX_VALUE, MatLib.abs(i));
    //noinspection NumericOverflow
    assertEquals(Integer.MAX_VALUE + 1, MatLib.abs(j));
    new Random().ints(1000).forEach(e -> assertEquals(Math.abs(e), MatLib.abs(e)));
  }

  @Test
  public void abs1() throws Exception {
    long i = Long.MAX_VALUE;
    long j = Long.MIN_VALUE;

    assertEquals(Long.MAX_VALUE, MatLib.abs(i));

    // Overflows because the min value has a larger
    // absolute value than the max value of a long
    //noinspection NumericOverflow
    assertEquals(Long.MAX_VALUE + 1, MatLib.abs(j));

    // Random doubles
    new Random().longs(1000).forEach(e -> assertEquals(Math.abs(e), MatLib.abs(e)));
  }

  @Test
  public void abs2() throws Exception {
    float i = Float.MAX_VALUE;
    float j = Float.MIN_VALUE;

    assertEquals(Float.MAX_VALUE, MatLib.abs(i), 0.0001);
    assertEquals(Float.MIN_VALUE, MatLib.abs(j), 0.0001);

    for(int l = 0; l < 1000; l++) {
      float n = new Random().nextFloat() * 2 - 1;
      assertEquals(Math.abs(n), MatLib.abs(n), 0.0001);
    }
  }

  @Test
  public void abs3() throws Exception {
    double i = Double.MAX_VALUE;
    double j = Double.MIN_VALUE;

    assertEquals(Double.MAX_VALUE, MatLib.abs(i), 0.0001);
    assertEquals(Double.MIN_VALUE, MatLib.abs(j), 0.0001);

    new Random().doubles(1000).forEach(e -> {
      e = 2 * e - 1;
      assertEquals(Math.abs(e), MatLib.abs(e), 0.0001);
    });
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
}