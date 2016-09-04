package com.adr.matlib.exception;

public class NonConformableMatrixException extends Exception {

    public NonConformableMatrixException() {
      this("Incorrect matrix sizes");
    }

    public NonConformableMatrixException(String msg) {
      super(msg);
    }
}

