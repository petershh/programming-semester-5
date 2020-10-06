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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matrixio.h"
#include "matrixlib.h"

int main(int argc, char **argv) {
    int n, m, k, result;
    double *matrix, *inverse;
    clock_t begin, end;
    int exit_code = 0;
    char *filename = NULL;
    if((argc < 4) || (argc > 5)) {
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
    if(sscanf(argv[3], "%d", &k) != 1) {
        exit_code = 1;
        goto final;
    }
    if(k < 0 || k > 4 || n < 1 || m < 1) {
        exit_code = 1;
        goto final;
    }
    if(argc == 5) {
        if(k != 0){
            exit_code = 1;
            goto final;
        }
        filename = argv[4];
    }

    matrix = (double*)malloc(n * n * sizeof(double));
    if(!matrix) {
        fprintf(stderr, "ERROR: not enough memory!");
        exit_code = 2;
        goto final;
    }
    inverse = (double*)malloc(n * n * sizeof(double));
    if(!inverse) {
        fprintf(stderr, "ERROR: not enough memory!");
        exit_code = 3;
        goto free_matrix;
    }

    if(read_matrix(matrix, n, k, filename)) {
        exit_code = 4;
        goto free_inverse;
    }

    printf("Original matrix:\n");
    print_matrix(matrix, n, n, m);
    printf("\n");

    begin = clock();
    result = invert_matrix(matrix, inverse, n);
    end = clock();

    if(result) {
        fprintf(stderr, "ERROR: matrix is not invertible\n");
        exit_code = 5;
        goto free_inverse;
    }

    printf("Inverted matrix:\n");
    print_matrix(inverse, n, n, m);
    printf("\n");

    if(read_matrix(matrix, n, k, filename)){
        exit_code = 4;
        goto free_inverse;
    }

    printf("Discrepancy: %e\n", discrepancy(matrix, inverse, n));
    printf("Time used to compute: %.2lf seconds\n", (double)(end - begin)
        / CLOCKS_PER_SEC);

    free_inverse:
    free(inverse);
    free_matrix:
    free(matrix);
    final:
    return exit_code;
}