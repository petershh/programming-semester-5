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

#include "matrixlib.h"
#include "common.h"

int invert_matrix(double *matrix, double *result, int order) {
	double s, norm1, norm2_square, tmp;

	// Generate the identity matrix

	memset(result, 0, (size_t)order * order * sizeof(double));
	for(int i = 0; i < order; i++)
		result[COORD(i, i, order)] = 1.0;
	
	// Cast the matrix to upper triangular type
	for(int i = 0; i < order; i++) {
		s = 0.0;
		for(int j = i + 1; j < order; j++){
			s += SQUARE(matrix[COORD(i, j, order)]);
		}

		norm1 = sqrt(SQUARE(matrix[COORD(i, i, order)]) + s);

		if(norm1 < EPS) {
			return 1; // non-invertible matrix
		}

		if(s < EPS) {
			continue; // nothing to do there
		}

		matrix[COORD(i, i, order)] -= norm1;
		norm2_square = SQUARE(matrix[COORD(i, i, order)]) + s;

		norm2_square = 2.0 / norm2_square;

		// Vector of reflection is ready, now we need to operate on matrices

		for(int j = i + 1; j < order; j++) {
			s = 0.0;
			for(int k = i; k < order; k++) {
				s += matrix[COORD(i, k, order)] * matrix[COORD(j, k, order)];
			}

			s *= norm2_square;
			for(int k = i; k < order; k++) {
				matrix[COORD(j, k, order)] -= s * matrix[COORD(i, k, order)];
			}
		}

		for(int j = 0; j < order; j++) {
			s = 0.0;
			for(int k = i; k < order; k++) {
				s += matrix[COORD(i, k, order)] * result[COORD(j, k, order)];
			}

			s *= norm2_square;
			for(int k = i; k < order; k++) {
				result[COORD(j, k, order)] -= s * matrix[COORD(i, k, order)];
			}
		}

		// Finalize: set the i-th subcolumn of matrix
		matrix[COORD(i, i, order)] = norm1;
	}

	// Back substitution of Gaussian method
	// We know that the matrix is inversible at the moment
	// Note: no action is required on matrix

	for(int i = order - 1; i >= 0; i--) {
		// Divide i-th row of result by matrix[i, i]

		s = matrix[COORD(i, i, order)];
		for(int j = 0; j < order; j++) {
			result[COORD(j, i, order)] /= s;
		}

		// Substract i-th row of result multiplied by matrix[j, i] from
		// j-th row of result for j = 0, ..., i - 1
		// But do it in the column-first order
		for(int j = 0; j < order; j++) {
			tmp = result[COORD(j, i, order)];
			for(int k = 0; k < i; k++) {
				result[COORD(j, k, order)] -= tmp * matrix[COORD(i, k, order)];
			}
		}
	}

	// And... here we go
	return 0;
}

double discrepancy(double *matrix, double *result, int order) {
	double product_elem = 0.0;
	double norm_square = 0.0;

	for(int i = 0; i < order; i++) {
		for(int j = 0; j < order; j++) {
			product_elem = 0.0;
			for(int k = 0; k < order; k++) {
				product_elem += matrix[COORD(k, i, order)] *
					result[COORD(j, k, order)];
			}

			norm_square += SQUARE(product_elem - (double)(i == j));
		}
	}

	return sqrt(norm_square);
}
