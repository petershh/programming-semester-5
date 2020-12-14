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

#include "common.h"
#include "matrixio.h"

int read_matrix(double *matrix, int order, int formula_number,
    char *filename) {
    if(filename) {
        FILE *fin = fopen(filename, "r");
        int result = 0;
        if(!fin) {
            perror("ERROR: failed to open file");
            return 1;
        }
        for(int i = 0; i < order; i++) {
            for(int j = 0; j < order; j++) {
                result = fscanf(fin, "%lf", matrix + COORD(j, i, order));
                if(result != 1) {
                    if(result == EOF) {
                        fprintf(stderr,
                            "ERROR: unexcepted EOF while reading matrix\n");
                    } else {
                        fprintf(stderr,
                            "ERROR: got invalid data while reading matrix\n");
                    }
                    fclose(fin);
                    return 1;
                }
            }
        }
        fclose(fin);
    } else {
        for(int i = 0; i < order; i++) {
            for(int j = 0; j < order; j++) {
                matrix[COORD(j, i, order)] = f(order, formula_number, i + 1,
                    j + 1);
            }
        }
    }
    return 0;
}

void print_matrix(double *matrix, int height, int width, int max_cols_rows) {
    int print_limit_x = MIN(width, max_cols_rows);
    int print_limit_y = MIN(height, max_cols_rows);
    for(int i = 0; i < print_limit_y; i++) {
        for(int j = 0; j < print_limit_x; j++) {
            printf(" %10.3e", matrix[COORD(j, i, height)]);
        }
        printf("\n");
    }
}