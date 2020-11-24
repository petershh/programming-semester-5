/*
 *  Copyright 2020 Peter Shkenev
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fenv.h>

#include "matrixio.h"
#include "matrixlib.h"

int main(int argc, char **argv) {
    int n, m, k;
    double *matrix, *eigenvalues, eps;
    clock_t begin, end;
    int exit_code = 0;
    char *filename = NULL;

//    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    if((argc < 5) || (argc > 6)) {
        exit_code = 1;
        goto final;
    }
    if(sscanf(argv[1], "%d", &n) != 1) {
        exit_code = 1;
        goto final;
    }
    if(sscanf(argv[2], "%d", &m) != 1) {
        exit_code = 1;
        goto final;
    }
    if(sscanf(argv[3], "%lf", &eps) != 1) {
        exit_code = 1;
        goto final;
    }
    if(sscanf(argv[4], "%d", &k) != 1) {
        exit_code = 1;
        goto final;
    }
    if(k < 0 || k > 4 || n < 1 || m < 1 || eps < 0.0) {
        exit_code = 1;
        goto final;
    }
    if(argc == 6) {
        if(k != 0){
            exit_code = 1;
            goto final;
        }
        filename = argv[5];
    }

    matrix = (double*)malloc(n * n * sizeof(double));
    if(!matrix) {
        fprintf(stderr, "ERROR: not enough memory!");
        exit_code = 2;
        goto final;
    }
    eigenvalues = (double*)malloc(n * sizeof(double));
    if(!eigenvalues) {
        fprintf(stderr, "ERROR: not enough memory!");
        exit_code = 3;
        goto free_matrix;
    }

    if(read_matrix(matrix, n, k, filename)) {
        exit_code = 4;
        goto free_eigenvalues;
    }

    printf("Original matrix:\n");
    print_matrix(matrix, n, n, m);
    printf("\n");

    begin = clock();
    get_eigenvalues(matrix, eigenvalues, n, eps);
    end = clock();

    printf("Eigenvalues:\n");
    print_matrix(eigenvalues, 1, n, n);
    printf("\n");

    if(read_matrix(matrix, n, k, filename)){
        exit_code = 4;
        goto free_eigenvalues;
    }

    printf("Residual 1: %e\n", residual1(matrix, eigenvalues, n));
    printf("Residual 2: %e\n", residual2(matrix, eigenvalues, n));
    printf("Time used to compute: %.2lf seconds\n", (double)(end - begin)
        / CLOCKS_PER_SEC);

    free_eigenvalues:
    free(eigenvalues);
    free_matrix:
    free(matrix);
    final:
    return exit_code;
}
