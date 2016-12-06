package com.adr.matlib;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Main {
    public static void main(String[] args) {
        List<String> pulseArrayList = readDoubleFile("pulse.txt");
        double[] pulseArray = new double[pulseArrayList.size()];

        for (int i = 0; i < pulseArray.length; i++) {
            pulseArray[i] = Double.valueOf(pulseArrayList.get(i));
        }

        List<String> receivingPulseArrayList = readDoubleFile("signal.txt");
        double[] receivingArray = new double[receivingPulseArrayList.size()];

        for (int i = 0; i < receivingArray.length; i++) {
            receivingArray[i] = Double.valueOf(receivingPulseArrayList.get(i));
        }

        /*for (double aReceivingArray : receivingArray) {
            System.out.println(aReceivingArray);
        }*/

        double[] correlation = MatLib.normalizedCrossCorrelatiton(pulseArray, receivingArray);

        for (double aCorrelation : correlation) {
            System.out.println(aCorrelation);
        }

        /*List<String> originalArrayList = readDoubleFile("original.txt");
        double[] originalArray = new double[originalArrayList.size()];
        for(int i = 0; i < originalArrayList.size(); i++) {
            originalArray[i] = Double.valueOf(originalArrayList.get(i));
        }

        double[] fftConvo = MatLib.fftConvolution(originalArray, 10);

        for (int i = 0; i < fftConvo.length; i++) {
            System.out.println(fftConvo[i]);
        }*/
    }

    private static List<String> readDoubleFile(String file) {
        String fileName = "resources/" + file;
        BufferedReader br;
        List<String> fileData = new ArrayList<>();

        try {
            br = new BufferedReader(new FileReader(new File(fileName)));
            String line;
            while ((line=br.readLine()) != null) {
                fileData.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return fileData;
    }
}
