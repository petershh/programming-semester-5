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

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "matrixlib.h"
#include "matrixio.h"
#include "common.h"

int get_eigenvalues(double *matrix, double *values, int order, double eps) {
    double temp1, temp2, temp3, temp4, temp5;
    double cos_phi = 0.0, sin_phi = 0.0;
    double cos_phi1 = 0.0, sin_phi1 = 0.0;
    double norm = infinity_norm(matrix, order);

    double *main_diag = matrix;
    double *lower_diag = matrix + order;

    memset(values, 0, order * sizeof(int));
    
    // Cast to three-diagonal type
    for(int i = 0; i < order - 2; i++) {
        // Working with vector (matrix[i, i + 1], ..., matrix[i, order])
        for(int j = i + 2; j < order; j++) {
            if(fabs(matrix[COORD(i, j, order)]) < 1e-16) {
                continue;
            }

            // Prepare matrix of rotation T = T(i + 1, j)
            temp1 = sqrt(SQUARE(matrix[COORD(i + 1, i, order)]) +
                SQUARE(matrix[COORD(i, j, order)]));
            temp2 = 1 / temp1;
            cos_phi = matrix[COORD(i + 1, i, order)] * temp2;
            sin_phi = -matrix[COORD(i, j, order)] * temp2;

            // We know what will happen with i-th column and i-th row
            // after multiplication
            matrix[COORD(i, i + 1, order)] = temp1;
            matrix[COORD(i + 1, i, order)] = temp1;
            matrix[COORD(i, j, order)] = 0.0;
            matrix[COORD(j, i, order)] = 0.0;


            // Multiply matrix by T from left
            for(int k = i + 1; k < order; k++) {
                // We need to preserve matrix[k, i + 1] value
                temp1 =
                    matrix[COORD(i + 1, k, order)] * cos_phi -
                    matrix[COORD(j, k, order)] * sin_phi;
                matrix[COORD(j, k, order)] = 
                    matrix[COORD(i + 1, k, order)] * sin_phi +
                    matrix[COORD(j, k, order)] * cos_phi;
                matrix[COORD(i + 1, k, order)] = temp1;
            }


            // "Multiply" matrix by T* from right
            // As stated in lemma 14.2.1 of The Book, we have only to
            // compute four elements, others can be relocated
            temp1 = 
                matrix[COORD(i + 1, i + 1, order)] * cos_phi-
                matrix[COORD(i + 1, j, order)] * sin_phi;
            temp2 = 
                matrix[COORD(i + 1, i + 1, order)] * sin_phi +
                matrix[COORD(i + 1, j, order)] * cos_phi;
            matrix[COORD(j, j, order)] = 
                matrix[COORD(j, i + 1, order)] * sin_phi +
                matrix[COORD(j, j, order)] * cos_phi;
            matrix[COORD(i + 1, i + 1, order)] = temp1;
            matrix[COORD(i + 1, j, order)] = temp2;
            matrix[COORD(j, i + 1, order)] = temp2;

			for(int k = i + 2; k < j; k++) {
				matrix[COORD(k, i + 1, order)] =
					matrix[COORD(i + 1, k, order)];
				matrix[COORD(k, j, order)] =
					matrix[COORD(j, k, order)];
			}

			for(int k = j + 1; k < order; k++) {
				matrix[COORD(k, i + 1, order)] =
					matrix[COORD(i + 1, k, order)];
				matrix[COORD(k, j, order)] =
					matrix[COORD(j, k, order)];
			}
        }
    }

    // relocate elements for easier code and speed
    for(int i = 1; i < order - 1; i++) {
        main_diag[i] = matrix[COORD(i, i, order)];
        lower_diag[i] = matrix[COORD(i + 1, i, order)];
    }
    main_diag[order - 1] = matrix[COORD(order - 1, order - 1, order)];
    lower_diag[order - 1] = 0.0;

    // Obtain eigenvalues
    for(int i = order - 1; i > 1; i--) {
        while(fabs(lower_diag[i - 1]) >= eps * norm) {
            temp1 = main_diag[i];  // Shift
            for(int j = 0; j <= i; j++) {
                main_diag[j] -= temp1;
            }
            //printf("\r");
            //printf("%lf", lower_diag[i - 1]);

            // QR decomposition (and computation of RQ at the same time)
			if(fabs(lower_diag[0]) >= 1e-16) { 
				temp2 = sqrt(SQUARE(main_diag[0]) +
					SQUARE(lower_diag[0]));
				temp3 = 1 / temp2;

				cos_phi = main_diag[0] * temp3;
				sin_phi = -lower_diag[0] * temp3;

				main_diag[0] = temp2;

				temp4 = lower_diag[0] * cos_phi -
					main_diag[1] * sin_phi; // 0-th elem of upper diag
				main_diag[1] = lower_diag[0] * sin_phi +
					main_diag[1] * cos_phi;

				lower_diag[0] = temp4; // Upper diag

            	temp5 = lower_diag[1] * cos_phi; // 1-th elem of upper diag

				cos_phi1 = cos_phi; sin_phi1 = sin_phi;
			} else {
				cos_phi1 = cos_phi = 1.0;
				sin_phi1 = sin_phi = 0.0;
				temp5 = lower_diag[1];
			}

            for(int k = 1; k <= i - 1; k++) {
				if(fabs(lower_diag[k]) >= 1e-16){
					temp2 = sqrt(SQUARE(main_diag[k]) +
						SQUARE(lower_diag[k]));
					temp3 = 1 / temp2;

					cos_phi = main_diag[k] * temp3;
					sin_phi = -lower_diag[k] * temp3;

					// Multiply by T(k, k + 1) from left
					main_diag[k] = temp2;
					temp4 = temp5 * cos_phi-
						main_diag[k + 1] * sin_phi;
					main_diag[k + 1] = temp5 * sin_phi +
						main_diag[k + 1] * cos_phi;
					lower_diag[k] = temp4;

					temp5 = lower_diag[k + 1] * cos_phi;
				} else {
					cos_phi = 1.0; sin_phi = 1.0;
					temp5 = lower_diag[k + 1];
				}

                // Multiply by T(k - 1, k)* from right
                main_diag[k - 1] = main_diag[k - 1] * cos_phi1 -
                    lower_diag[k - 1] * sin_phi1;

                // Corresponding upper diagonal element will be just
                // lower_diag[k - 1] * cos_phi
                lower_diag[k - 1] = - main_diag[k] * sin_phi1;
                main_diag[k] *= cos_phi1;

                cos_phi1 = cos_phi;
                sin_phi1 = sin_phi;
            }

            // Last iteration
            main_diag[i - 1] = main_diag[i - 1] * cos_phi1 - 
                lower_diag[i - 1] * sin_phi1;
            
            lower_diag[i - 1] = -main_diag[i] * sin_phi1;
            main_diag[i] = main_diag[i] * cos_phi1;

            // Shift back
            for(int j = 0; j <= i; j++) {
                main_diag[j] += temp1;
            }
        }
        values[i] = main_diag[i];
		lower_diag[i - 1] = 0.0;
    }
    // Last two eigenvalues are found as solutions of quadratic equation
    temp1 = sqrt(SQUARE(main_diag[0] - main_diag[1])
        + 4 * SQUARE(lower_diag[0]));
    values[1] = (main_diag[0] + main_diag[1] - temp1) / 2.0;
    values[0] = (main_diag[0] + main_diag[1] + temp1) / 2.0;
    return 0;
}

double residual1(double *matrix, double *eigenvalues, int order) {
    double result = 0.0;
    for(int i = 0; i < order; i++) {
        result += (matrix[COORD(i, i, order)] - eigenvalues[i]);
    }
    return fabs(result);
}

double residual2(double *matrix, double *eigenvalues, int order) {
    double result = 0.0, result1 = 0.0;
    for(int i = 0; i < SQUARE(order); i++) {
        result += SQUARE(matrix[i]);
    }
	result = sqrt(result);
    for(int i = 0; i < order; i++) {
        result1 += SQUARE(eigenvalues[i]);
    }
	result1 = sqrt(result1);
    return fabs(result - result1);
}

double infinity_norm(const double *matrix, int order) {
    double t;
    double result = -1.0;

    for(int i = 0; i < order; i++) {
        t = 0.0;
        for(int j = 0; j < order; j++) {
            t += fabs(matrix[COORD(i, j, order)]);
        }
        result = MAX(result, t);
    }

    return result;
}

