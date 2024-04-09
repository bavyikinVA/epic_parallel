#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>

double** morlet_wavelet(double** data, int data_rows, int data_cols, double* weight, int weight_size) 
{    
    //double** coef = (double**)malloc(weight_size * sizeof(double*));
    double*** result_matrix = new double** [weight_size * data_cols * data_cols];
    for (int i = 0; i < weight_size; i++) 
    {
        double** coef = new double*[weight_size];
        for (int l = 0; l < data_rows; l++) 
        {
            coef[i] = new double*[data_size];
            double* coefRow = new double[data_size];
            //omp_set_num_threads(4);
            //#pragma omp parallel for schedule(static)
            for (int j = 0; j < data_size; j++) 
            {
                double w0 = 0;
                for (int k = 1; k < data_size; k++) {
                    double t = (double)(k - j) / weight[i];
                    w0 = w0 + data[k] * 0.75 * exp(-(t * t) / 2) * cos(6 * t);
                }
                coefRow[j] = w0 / sqrt(weight[i]);

            }
        coef[i] = coefRow;
    }

    return coef;
}

int main() {
    std::string input_path;
    std::cout << "Enter file path:" << std::endl;

    std::cin >> input_path;
    std::cout << std::endl;

    for (size_t i = 0; i < input_path.length(); ++i) 
    {
        if (input_path[i] == '/' || input_path[i] == '\\') {
            input_path.replace(i, 1, "\\\\");
            ++i;
        }
    }

    std::ifstream fin(input_path);
    if (!fin.is_open()) 
    {
        std::cerr << "Error opening file: " << input_path << std::endl;
        return 1;
    }

    int rows = 0;
    int cols = 0;

    // Read matrix dimensions
    fin >> rows >> cols;

    double** matrix = new double* [rows];
    for (int i = 0; i < rows; ++i) {
        matrix[i] = new double[cols];
    }

    // Read matrix elements
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            fin >> matrix[i][j];
        }
    }
    fin.close();

    std::cout << "select an option: 1 - list of wavelet transform scales, 2 - interval of the wavelet transform scale:" << std::endl;
    int selector;
    int iter = 0;
    double* wavelet_scales = new double[1];
    std::cin >> selector;
    if (selector == 1) 
    {
        std::cout << "Enter the number of items:" << std::endl;
        std::cin >> iter;
        double* scales2 = new double[iter];
        for (int i = 0; i < iter; i++)
        {
            std::cin >> scales2[i];

        }
        wavelet_scales = new double[iter];
        for (int i = 0; i < iter; i++) 
        {
            wavelet_scales[i] = scales2[i];
        }
    }



    /*double widths[] = {}
    for (int i = 0; i < WEIGHT_SIZE; ++i) {
        widths[i] = i + 1;
    }*/

    double time = (double)clock() / CLOCKS_PER_SEC;
    printf("start_morlet");
    double** result = morlet_wavelet(matrix, rows, cols, wavelets_scales, iter);
    double time_diff = (((double)clock()) / CLOCKS_PER_SEC) - time;

    std::ofstream fout("C:\\Users\\bavyk\\source\\repos\\epic_parallel\\epic_parallel\\test.txt");

    // Вывод результатов
    double time2 = ((double)clock()) / CLOCKS_PER_SEC;
    for (int i = 0; i < WEIGHT_SIZE; i++)
    {
        for (int j = 0; j < DATA_SIZE; j++)
        {
            fout << std::to_string(result[i][j]) << " ";
        }

        fout << std::endl;
    
    }
    fout.close();
    double time_diff2 = (((double)clock()) / CLOCKS_PER_SEC) - time2;
    printf("The elapsed time is %lf seconds\n", time_diff);
    printf("The elapsed time is %lf seconds\n", time_diff2);

    // Освобождение памяти
    for (int i = 0; i < WEIGHT_SIZE; i++) 
    {
        free(result[i]);
    }
    free(result);
    system("pause");
    return 0;
}