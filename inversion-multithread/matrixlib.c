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
#include <pthread.h>
#include <assert.h>

#include "matrixlib.h"
#include "common.h"

int invert_matrix(double *matrix, double *result, int order, int thread_id,
					int threads_amount) {
	double s, norm1, norm2, tmp;
	int work_range_start, work_range_end;

	synchronize(threads_amount);

	// Generate the identity matrix
	
	if(thread_id == 0) {
		memset(result, 0, (size_t)order * order * sizeof(double));
		for(int i = 0; i < order; i++)
			result[COORD(i, i, order)] = 1.0;
	}
	
	synchronize(threads_amount);
	// Cast the matrix to upper triangular type
	for(int i = 0; i < order; i++) {
		s = 0.0;
		for(int j = i + 1; j < order; j++) {
			s += SQUARE(matrix[COORD(i, j, order)]);
		}

		norm1 = sqrt(SQUARE(matrix[COORD(i, i, order)]) + s);

		if(norm1 < EPS) {
			synchronize(threads_amount);	
			return 1; // non-invertible matrix
		}

		if(s < EPS) {
			synchronize(threads_amount);
			continue; // nothing to do there
		}

		if(thread_id == 0) {
			matrix[COORD(i, i, order)] -= norm1;
			norm2 = sqrt(SQUARE(matrix[COORD(i, i, order)]) + s);

			norm2 = 1.0 / norm2;
			for(int j = i; j < order; j++) {
				matrix[COORD(i, j, order)] *= norm2;
			}
		}

		// Vector of reflection is ready, now we need to operate on matrices
		synchronize(threads_amount);

		work_range_start = ((order - i - 1) * thread_id) / threads_amount +
			i + 1;
		work_range_end = ((order - i - 1) * (thread_id + 1)) / threads_amount +
		   	i + 1;

		for(int j = work_range_start; j < work_range_end; j++) {
			s = 0.0;
			for(int k = i; k < order; k++) {
				s += matrix[COORD(i, k, order)] * matrix[COORD(j, k, order)];
			}

			s *= 2.0;
			for(int k = i; k < order; k++) {
				matrix[COORD(j, k, order)] -= s * matrix[COORD(i, k, order)];
			}
		}

		work_range_start = (order * thread_id) / threads_amount;
		work_range_end = (order * (thread_id + 1)) / threads_amount;

		for(int j = work_range_start; j < work_range_end; j++) {
			s = 0.0;
			for(int k = i; k < order; k++) {
				s += matrix[COORD(i, k, order)] * result[COORD(j, k, order)];
			}

			s *= 2.0;
			for(int k = i; k < order; k++) {
				result[COORD(j, k, order)] -= s * matrix[COORD(i, k, order)];
			}
		}

		synchronize(threads_amount);

		// Finalize: set the i-th subcolumn of matrix
		if(thread_id == 0) {
			matrix[COORD(i, i, order)] = norm1;
		}
	}

	// Back substitution of Gaussian method
	// We know that the matrix is inversible at the moment
	// Note: no action is required on matrix

	for(int i = order - 1; i >= 0; i--) {
		// Divide i-th row of result by matrix[i, i]

		s = 1 / matrix[COORD(i, i, order)];
		if(thread_id == 1) {
			for(int j = 0; j < order; j++) {
				result[COORD(j, i, order)] *= s;
			}
		}
		synchronize(threads_amount);

		// Substract i-th row of result multiplied by matrix[j, i] from
		// j-th row of result for j = 0, ..., i - 1
		// But do it in the column-first order

		work_range_start = (order * thread_id) / threads_amount;
		work_range_end = (order * (thread_id + 1)) / threads_amount;

		for(int j = work_range_start; j < work_range_end; j++) {
			tmp = result[COORD(j, i, order)];
			for(int k = 0; k < i; k++) {
				result[COORD(j, k, order)] -= tmp * matrix[COORD(i, k, order)];
			}
		}
		synchronize(threads_amount);
	}

	// And... here we go
	return 0;
}

double residual(double *matrix, double *result, int order, int thread_id,
				int threads_amount) {
	double product_elem = 0.0;
	double norm_square = 0.0;
	int work_range_start = (order * thread_id) / threads_amount;
	int work_range_end = (order * (thread_id + 1)) / threads_amount;

	for(int i = work_range_start; i < work_range_end; i++) {
		for(int j = 0; j < order; j++) {
			product_elem = 0.0;
			for(int k = 0; k < order; k++) {
				product_elem += matrix[COORD(k, i, order)] *
					result[COORD(j, k, order)];
			}

			norm_square += SQUARE(product_elem - (double)(i == j));
		}
	}

	return norm_square;
}
