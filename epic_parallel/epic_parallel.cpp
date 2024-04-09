#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>

#define DATA_SIZE 1001
#define WEIGHT_SIZE 100

double** morlet_wavelet(double* data, int data_size, double* weight, int weight_size) 
{    
    double** coef = (double**)malloc(weight_size * sizeof(double*));
    for (int i = 0; i < weight_size; i++) {
        coef[i] = (double*)malloc(data_size * sizeof(double));
        double* coefRow = (double*)malloc(data_size * sizeof(double));
        omp_set_num_threads(4);
        #pragma omp parallel for schedule(static)
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
    std::ifstream fin("C:\\Users\\bavyk\\source\\repos\\epic_parallel\\epic_parallel\\delta_.txt");
    int size = 1001;
    double delta[1001] = {};
    //delta[501] = {1.0};
    // Reading the array elements from the file 
    for (int i = 0; i < size; i++)
    {
        fin >> delta[i];
    }
    fin.close();
    double widths[WEIGHT_SIZE] = {};
    for (int i = 0; i < WEIGHT_SIZE; ++i) {
        widths[i] = i + 1;
    }

    double time = (double)clock() / CLOCKS_PER_SEC;
    printf("start_morlet");
    double** result = morlet_wavelet(delta, DATA_SIZE, widths, WEIGHT_SIZE);
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