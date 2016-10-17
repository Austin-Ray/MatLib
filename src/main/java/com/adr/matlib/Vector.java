package com.adr.matlib;

public class Vector {

  double x;
  double y;

  public Vector(double x, double y) {
    this.x = x;
    this.y = y;
  }

  public double getX() {
    return x;
  }

  public double getY() {
    return y;
  }

  public double[][] matrix() {
    return new double[][] {{x}, {y}};
  }
}
