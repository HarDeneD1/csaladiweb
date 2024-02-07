package linear.algebra;

import java.util.Arrays;

public class GaussianElimination {

    private final int numRows;
    private final int numCols;
    private final double[][] matrix;

    public GaussianElimination(int numRows, int numCols, double[][] matrix) {
       
        this.numRows = numRows;
        this.numCols = numCols;
        this.matrix = new double[numRows][numCols];
       
        for (int i = 0; i < numRows; i++) {
            this.matrix[i] = Arrays.copyOf(matrix[i], numCols);
        }
    }

    public int getNumRows() {
        return numRows;
    }

    public int getNumCols() {
        return numCols;
    }

    public double[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(double[][] newMatrix) {
        if (newMatrix.length != numRows || newMatrix[0].length != numCols) {
            throw new IllegalArgumentException("New matrix must have the same number of rows and columns as the original matrix.");
        }
        for (int i = 0; i < numRows; i++) {
            matrix[i] = Arrays.copyOf(newMatrix[i], numCols);
        }
    }

    public void rowEchelonForm() {
        int lead = 0;
        for (int r = 0; r < numRows; r++) {
            if (numCols <= lead) {
                return;
            }
            int i = r;
            while (matrix[i][lead] == 0) {
                i++;
                if (numRows == i) {
                    i = r;
                    lead++;
                    if (numCols == lead) {
                        return;
                    }
                }
            }
            double[] temp = matrix[i];
            matrix[i] = matrix[r];
            matrix[r] = temp;

            double lv = matrix[r][lead];
            for (int j = 0; j < numCols; j++) {
                matrix[r][j] /= lv;
            }
            for (int k = 0; k < numRows; k++) {
                if (k != r) {
                    double mult = matrix[k][lead];
                    for (int j = 0; j < numCols; j++) {
                        matrix[k][j] -= mult * matrix[r][j];
                    }
                }
            }
            lead++;
        }
    }
}
