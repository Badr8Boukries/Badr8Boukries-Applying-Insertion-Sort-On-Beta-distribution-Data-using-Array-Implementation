package org.example;



import org.apache.commons.math3.distribution.BetaDistribution;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import java.text.DecimalFormat;

import javax.swing.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

public class BetaDistributionMatrices {

    // Function to generate an array of random numbers from the beta distribution
    static double[] betaDataGenerate(double mean, double variance, int dataSize) {
        if (dataSize <= 0) {
            throw new IllegalArgumentException("Data size must be Positif");
        }
        double[] data = new double[dataSize];
        double alpha = mean * (mean * (1 - mean) / variance - 1);
        double beta = (1 - mean) * (mean * (1 - mean) / variance - 1);
        BetaDistribution betaDistribution = new BetaDistribution(alpha, beta);
        for (int i = 0; i < dataSize; i++) {
            data[i] = betaDistribution.sample();
        }
        return data;
    }

    // Function to sort the data using the insertion sort algorithm
    public static void insertionSort(double[] arr) {
        for (int i = 1; i < arr.length; i++) {
            double key = arr[i];
            int j = i - 1;
            while (j >= 0 && arr[j] > key) {
                arr[j + 1] = arr[j];
                j--;
            }
            arr[j + 1] = key;
        }
    }

    // Fonction pour chauffer le JVM
    public static void warmUpJVM() {
        System.out.println("Chauffage du JVM en cours...");
        for (int i = 0; i < 100; i++) {
            double[] randomNumbers = betaDataGenerate(0.5, 100, 1000);
            insertionSort(randomNumbers);
        }
        System.out.println("Chauffage du JVM termine.");
    }

    // Function to warm up the JVM
    public static double[][][] createBetaMatrices(double[] means, double[] variances, int initialDataSize) {
        int numMeans = means.length;
        int numVariances = variances.length;
        double[][][] matrices = new double[numMeans][numVariances][];

        for (int i = 0; i < numMeans; i++) {
            for (int j = 0; j < numVariances; j++) {
                double[] data = betaDataGenerate(means[i], variances[j], initialDataSize);

                matrices[i][j] = data;

                // Store the sorted data in a CSV file for the average case
                storeDataInCSV(data, means[i], variances[j], initialDataSize, "Average case");


                double[] bestCaseData = Arrays.copyOf(data, data.length);
                insertionSort(bestCaseData);

                double[] worstCaseData = Arrays.copyOf(data, data.length);
                insertionSort(worstCaseData);
                reverseArray(worstCaseData);
            }
        }

        return matrices;
    }


    // Function to store the sorted data in a CSV file for the average case
    public static void storeDataInCSV(double[] data, double mean, double variance, int dataSize, String caseType) {
        // Crée le répertoire pour stocker les fichiers CSV s'il n'existe pas encore
        Path storedDataPath = Paths.get("storedData");
        if (!Files.exists(storedDataPath)) {
            try {
                Files.createDirectories(storedDataPath);
            } catch (IOException e) {
                System.err.println("Erreur : " + e.getMessage());
                return;
            }
        }

        // Crée le nom de fichier basé sur les paramètres fournis
        String fileName = storedDataPath + "/" + caseType.toLowerCase().replace(" ", "_") + "_case_results_datasize_" + dataSize + "_mean_" + mean + "_variance_" + variance + ".csv";

        // Copie les données et mesure le temps d'exécution pour les trier
        double[] sortedData = Arrays.copyOf(data, data.length);

        double executionTime = measureExecutionTime( sortedData);

        // Crée un objet DecimalFormat pour limiter la précision à 4 chiffres après la virgule
        DecimalFormat df = new DecimalFormat("#.####");

        try (FileWriter writer = new FileWriter(fileName)) {
            // Écrit l'entête du fichier CSV
            writer.append("mean       ,        variance       ,          dataSize    ,      executionTime     ,       sortedData\n");

            // Écrit les données formatées dans le fichier CSV
            writer.append(df.format(mean)).append("          ,     ");
            writer.append(df.format(variance)).append("        ,   ");
            writer.append(String.valueOf(dataSize)).append("          ,      ");
            writer.append(String.valueOf(executionTime)).append("   ms      , ");

            // Écrit les données triées formatées dans le fichier CSV
            for (int i = 0; i < sortedData.length; i++) {
                writer.append(df.format(sortedData[i]));
                if (i < sortedData.length - 1) {
                    writer.append("  ,  ");
                }
            }

            writer.append("\n");
        } catch (IOException e) {
            System.err.println("Erreur : " + e.getMessage());
        }
    }

    public static double calculateAverageExecutionTime(double[] data) {
        long totalExecutionTime = 0;

        // Effectuer le tri 10 fois et calculer le temps d'exécution total
        for (int i = 0; i < 2; i++) {
            // Créer une copie de l'array de données
            double[] copy = Arrays.copyOf(data, data.length);


            long startTime = System.nanoTime();
            insertionSort(copy);
            long endTime = System.nanoTime();


            totalExecutionTime += (endTime - startTime);
        }


        double averageExecutionTime = (totalExecutionTime / 2) / 1_000_000;

        return averageExecutionTime;
    }



    // Function to measure execution times for each combination of mean and variance in the 6 cases
    public static void measureExecutionTimes(double[][][] betaMatrices1, double[][][] betaMatrices2, double[][][] betaMatrices3,
                                             double[][][] betaMatrices4, double[][][] betaMatrices5, double[][][] betaMatrices6,
                                             double[] means, double[] variances, String caseType) {
        // Array of matrices for each data size
        double[][][][] betaMatricesList = {betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6};

        // Array to store the average execution times for each data size
        double[] meanExecutionTimesPerDataSize = new double[betaMatricesList.length];

        System.out.println("\nTemps d'exécution moyen pour le cas " + caseType + " :");

        // Array to store execution times
        double[] meanExecutionTimes = new double[means.length * variances.length];

        for (int i = 0; i < means.length; i++) {
            for (int j = 0; j < variances.length; j++) {
                double totalExecutionTime = 0;
                int count = 0;

                // Iterate through each data size and calculate the execution time for each combination of mean and variance
                for (int k = 0; k < betaMatricesList.length; k++) {
                    double[] data = Arrays.copyOf(betaMatricesList[k][i][j], betaMatricesList[k][i][j].length);

                    // If the case is "Average case", no transformation is necessary
                    // If the case is "Best case", sort the data
                    if (caseType.equals("Best case")) {
                        insertionSort(data);
                        perturbBestCaseData(data, 4000);
                    } else if (caseType.equals("Worst case")) {
                        // If the case is "Worst case", sort the data and reverse it
                        insertionSort(data);
                        reverseArray(data);
                    }

                    // Measure the execution time in double precision
                    double averageExecutionTime = calculateAverageExecutionTime(data);
                    totalExecutionTime += averageExecutionTime;
                    count++;

                    // Accumulate execution time for each data size
                    meanExecutionTimesPerDataSize[k] += averageExecutionTime;
                }

                // Calculate the average execution time for each combination of mean and variance
                double meanExecutionTime = totalExecutionTime / count;

                // Display the average execution time
                System.out.println("Mean: " + means[i] + ", Variance: " + variances[j] + ", Mean Execution Time: " + meanExecutionTime + " ms");

                // Store the average execution time in the array
                meanExecutionTimes[i * variances.length + j] = meanExecutionTime;
            }
        }

        // Calculate the average execution times for each data size
        System.out.println("\nAverage execution times for each data size (" + caseType + " case):");
        for (int k = 0; k < betaMatricesList.length; k++) {

            // Calculate the average execution time for each data size
            double meanExecutionTimePerDataSize = meanExecutionTimesPerDataSize[k] / (means.length * variances.length);
            System.out.println("Data size nk = " + (k + 1) + "0000, Average execution time: " + meanExecutionTimePerDataSize + " ms");
        }

        // Iterate through each data size and display the execution times for each combination of mean and variance
        if (caseType.equals("Average case")) {
            System.out.println("\nDetails of execution times for each combination of mean and variance for each data size (" + caseType + " case):");
            for (int nkIndex = 0; nkIndex < betaMatricesList.length; nkIndex++) {

                System.out.println("\nData size nk = " + (nkIndex + 1) + "0000:");
                for (int i = 0; i < means.length; i++) {
                    for (int j = 0; j < variances.length; j++) {

                        // Copy the data from the beta matrix for the current data size, specific mean, and variance
                        double[] data = Arrays.copyOf(betaMatricesList[nkIndex][i][j], betaMatricesList[nkIndex][i][j].length);

                        // Calculate the execution time in double precision for the current combination
                        double executionTime = measureExecutionTime(data);

                        // Display the details for each combination of mean and variance
                        System.out.println("Mean: " + means[i] + ", Variance: " + variances[j] + ", Execution time: " + executionTime + " ms");
                    }
                }
            }
        }
    }

    // Fonction pour visualiser l'impact de la variance sur le temps d'exécution avec un mean fixe
    public static void plotMeanImpact(double[][][] betaMatrices1, double[][][] betaMatrices2, double[][][] betaMatrices3,
                                      double[][][] betaMatrices4, double[][][] betaMatrices5, double[][][] betaMatrices6,
                                      double[] means, double[] variances, double varianceToUse, String[] cases) {
        double[][][][] betaMatricesList = {betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6};
        int[] dataSizes = {10000, 20000, 30000, 40000, 50000, 60000};
        int varianceIndex = findIndex(variances, varianceToUse);

        if (varianceIndex == -1) {
            System.out.println("Error: Variance not found.");
            return;
        }

        XYSeriesCollection dataset = new XYSeriesCollection();

        for (String caseType : cases) {
            XYSeries series = new XYSeries(caseType);

            for (int i = 0; i < means.length; i++) {
                double mean = means[i];
                long totalExecutionTime = 0;

                for (int j = 0; j < betaMatricesList.length; j++) {
                    double[] cell = Arrays.copyOf(betaMatricesList[j][i][varianceIndex], betaMatricesList[j][i][varianceIndex].length);

                    if (caseType.equals("Best case")) {
                        insertionSort(cell);
                        perturbBestCaseData(cell,4000);
                    } else if (caseType.equals("Worst case")) {
                        insertionSort(cell);
                        reverseArray(cell);
                    }

                    totalExecutionTime += (long) calculateAverageExecutionTime(cell);
                }

                double averageExecutionTime = (double) totalExecutionTime / betaMatricesList.length;
                series.add(mean, averageExecutionTime);
            }

            dataset.addSeries(series);
        }

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Impact of Mean on Execution Time",
                "Mean",
                "Execution Time (ms)",
                dataset
        );

        ChartPanel chartPanel = new ChartPanel(chart);
        JFrame frame = new JFrame("Impact of Mean on Execution Time");
        frame.setContentPane(chartPanel);
        frame.setSize(800, 600);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }


    // Fonction pour visualiser l'impact de la mean sur le temps d'exécution avec une variance fixe
    public static void plotVarianceImpact(double[][][] betaMatrices1, double[][][] betaMatrices2, double[][][] betaMatrices3,
                                          double[][][] betaMatrices4, double[][][] betaMatrices5, double[][][] betaMatrices6,
                                          double[] means, double[] variances, double meanToUse, String[] cases) {
        double[][][][] betaMatricesList = {betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6};
        int[] dataSizes = {10000, 20000, 30000, 40000, 50000, 60000};
        int meanIndex = findIndex(means, meanToUse);

        if (meanIndex == -1) {
            System.out.println("Error: Mean not found.");
            return;
        }

        XYSeriesCollection dataset = new XYSeriesCollection();

        for (String caseType : cases) {
            XYSeries series = new XYSeries(caseType);

            for (int i = 0; i < variances.length; i++) {
                double variance = variances[i];
                long totalExecutionTime = 0;

                for (int j = 0; j < betaMatricesList.length; j++) {
                    double[] cell = Arrays.copyOf(betaMatricesList[j][meanIndex][i], betaMatricesList[j][meanIndex][i].length);

                    if (caseType.equals("Best case")) {
                        insertionSort(cell);
                        perturbBestCaseData(cell,4000);
                    } else if (caseType.equals("Worst case")) {
                        insertionSort(cell);
                        reverseArray(cell);
                    }

                    totalExecutionTime += (long) calculateAverageExecutionTime(cell);
                }

                double averageExecutionTime = (double) totalExecutionTime / betaMatricesList.length;
                series.add(variance, averageExecutionTime);
            }

            dataset.addSeries(series);
        }

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Impact of Variance on Execution Time",
                "Variance",
                "Execution Time (ms)",
                dataset
        );

        ChartPanel chartPanel = new ChartPanel(chart);
        JFrame frame = new JFrame("Impact of Variance on Execution Time");
        frame.setContentPane(chartPanel);
        frame.setSize(800, 600);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }



    // Function to disturb the sorted data in the best case
    public static void perturbBestCaseData(double[] data, int numPerturbations) {
        int length = data.length;
        numPerturbations = Math.min(numPerturbations, length / 10);


        for (int i = 0; i < numPerturbations; i++) {
            // Choose two random indices to swap the values
            int index1 = (int) (Math.random() * length);
            int index2 = (int) (Math.random() * length);


            double temp = data[index1];
            data[index1] = data[index2];
            data[index2] = temp;
        }
    }




    // Function to reverse an array
    public static void reverseArray(double[] arr) {
        int start = 0;
        int end = arr.length - 1;
        while (start < end) {
            double temp = arr[start];
            arr[start] = arr[end];
            arr[end] = temp;
            start++;
            end--;
        }
    }

    // Function to measure the execution time of an insertion sort on a cell
    public static long measureExecutionTime(double[] array) {
        long startTime = System.nanoTime();
        insertionSort(array);
        long endTime = System.nanoTime();

        return (long) (((endTime - startTime) / 1_000_000));
    }



    public static void plotExecutionTime(double[][][] betaMatrices1, double[][][] betaMatrices2, double[][][] betaMatrices3,
                                         double[][][] betaMatrices4, double[][][] betaMatrices5, double[][][] betaMatrices6,
                                         double mean, double variance, double[] means, double[] variances, String[] cases) {
        // Array of matrices for each data size
        double[][][][] betaMatricesList = {betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6};
        int[] dataSizes = {1000, 2000, 30000, 40000, 50000, 60000};

        // Find the indices of mean and variance
        int meanIndex = findIndex(means, mean);
        int varianceIndex = findIndex(variances, variance);

        // Check if the indices are valid
        if (meanIndex == -1 || varianceIndex == -1) {
            System.out.println("Error: mean or variance not found.");
            return;
        }

        // Create a collection of data series
        XYSeriesCollection dataset = new XYSeriesCollection();

        // Iterate through each case (Average case, Best case, Worst case)
        for (int caseIndex = 0; caseIndex < cases.length; caseIndex++) {
            /// Create a new series for the current case
            XYSeries series = new XYSeries(cases[caseIndex]);

            /// Iterate through each data size
            for (int i = 0; i < dataSizes.length; i++) {

                /// Retrieve data from the beta matrix for the specific combination of mean and variance
                double[] cell = betaMatricesList[i][meanIndex][varianceIndex];

                // Appliquer les transformations de cas (average, best, worst)
                if (cases[caseIndex].equals("Best case")) {
                    insertionSort(cell);
                    perturbBestCaseData(cell,4000);
                } else if (cases[caseIndex].equals("Worst case")) {
                    insertionSort(cell);
                    reverseArray(cell);
                }

                // Mesure the execution time for the modified data
                long executionTime = (long) calculateAverageExecutionTime(cell);

                // Add the data size and execution time to the series
                series.add(dataSizes[i], executionTime);
            }


            // Add the series to the dataset
            dataset.addSeries(series);
        }

        JFreeChart chart = ChartFactory.createXYLineChart(
                "Execution Time for Different Data Sizes and Cases",
                "Data Size",
                "Execution Time (ms)",
                dataset
        );

        ChartPanel chartPanel = new ChartPanel(chart);
        JFrame frame = new JFrame("Execution Time Chart");
        frame.setContentPane(chartPanel);
        frame.setSize(800, 600);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }



    // Function to find the index of a value in an array
    public static int findIndex(double[] array, double value) {
        for (int i = 0; i < array.length; i++) {
            if (Double.compare(array[i], value) == 0) {
                return i;
            }
        }
        return -1;
    }

    public static void main(String[] args) {
        // Jvm Warm up
        warmUpJVM();

        double[] means = {0.1, 0.3, 0.5, 0.7, 0.9,0.99, 0.8, 0.6, 0.4, 0.2};
        double[] variances = {5, 10, 1500, 2000, 2500, 3000, 3500, 4000, 45000, 500000};
        int initialDataSize = 1000;

        String[] cases = {"Average case", "Best case", "Worst case"};

        // Create matrices for each combination of mean and variance
        double[][][] betaMatrices1 = createBetaMatrices(means, variances, initialDataSize);
        double[][][] betaMatrices2 = createBetaMatrices(means, variances,initialDataSize*2);
        double[][][] betaMatrices3 = createBetaMatrices(means, variances, initialDataSize*4);
        double[][][] betaMatrices4 = createBetaMatrices(means, variances, initialDataSize*8);
        double[][][] betaMatrices5 = createBetaMatrices(means, variances, initialDataSize*16);
        double[][][] betaMatrices6 = createBetaMatrices(means, variances, initialDataSize*32);

        // Measure the execution times
        measureExecutionTimes(betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6,
                means, variances, "Average case");
        measureExecutionTimes(betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6,
                means, variances, "Best case");
        measureExecutionTimes(betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6,
                means, variances, "Worst case");

        plotExecutionTime(betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6,
                0.5, 1500, means, variances, new String[]{"Average case", "Best case", "Worst case"});

        // Pour afficher l'impact de variance sur le temps d'exécution, avec un mean fixé à 0.5
        plotVarianceImpact(betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6,
                means, variances, 0.9, new String[]{"Average case", "Best case", "Worst case"});

// Pour afficher l'impact de mean sur le temps d'exécution, avec une variance fixée à 90
        plotMeanImpact(betaMatrices1, betaMatrices2, betaMatrices3, betaMatrices4, betaMatrices5, betaMatrices6,
                means, variances, 45000, new String[]{"Average case", "Best case", "Worst case"});




    }
}
