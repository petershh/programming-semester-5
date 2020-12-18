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
#include <math.h>
#include <pthread.h>

#include "common.h"
#include "matrixio.h"
#include "matrixlib.h"


pthread_mutex_t total_time_mutex = PTHREAD_MUTEX_INITIALIZER;
long int thread_total_time = 0L;
int inversion_result = 0;

pthread_mutex_t read_matrix_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t read_matrix_condvar = PTHREAD_COND_INITIALIZER;

struct thread_args {
	int thread_id;
	int threads_amount;
	double *matrix;
	double *inverse_matrix;
	int order;
	double residual_part;
};

void *thread_execute(void *p_args);

int main(int argc, char **argv) {
    int n, m, k, threads_amount;
    double *matrix, *inverse, residual_value;
    int exit_code = 0;
    char *filename = NULL;
	struct thread_args *args;
	pthread_t *threads;
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
    if(sscanf(argv[3], "%d", &k) != 1) {
        exit_code = 1;
        goto final;
    }
	if(sscanf(argv[4], "%d", &threads_amount) != 1) {
		exit_code = 1;
		goto final;
	}
    if(k < 0 || k > 4 || n < 1 || m < 1 || threads_amount < 1) {
        exit_code = 1;
        goto final;
    }
    if(argc == 6) {
        if(k != 0){
            exit_code = 1;
            goto final;
        }
        filename = argv[4];
    }

    matrix = (double*)malloc(n * n * sizeof(double));
    if(!matrix) {
        fprintf(stderr, "ERROR: not enough memory!\n");
        exit_code = 2;
        goto final;
    }

    inverse = (double*)malloc(n * n * sizeof(double));
    if(!inverse) {
        fprintf(stderr, "ERROR: not enough memory!\n");
        exit_code = 3;
        goto free_matrix;
    }

	args = (struct thread_args*)malloc(threads_amount *
			sizeof(struct thread_args));
	if(!args) {
		fprintf(stderr, "ERROR: not enough memory!\n");
		exit_code = 4;
		goto free_inverse;
	}

	threads = (pthread_t*)malloc(threads_amount * sizeof(pthread_t));
	if(!threads) {
		fprintf(stderr, "ERROR: not enough memory!\n");
		exit_code = 5;
		goto free_args;
	}

    if(read_matrix(matrix, n, k, filename)) {
        exit_code = 5;
        goto free_threads;
    }

	for(int i = 0; i < threads_amount; i++) {
		args[i].inverse_matrix = inverse;
		args[i].matrix = matrix;
		args[i].order = n;
		args[i].thread_id = i;
		args[i].threads_amount = threads_amount;
	}

    printf("Original matrix:\n");
    print_matrix(matrix, n, n, m);
    printf("\n");

	for(int i = 0; i < threads_amount; i++) {
		if(pthread_create(threads + i, NULL, thread_execute, args + i)) {
			fprintf(stderr, "ERROR: Cannot create threads!\n");
			exit_code = 7;
			goto free_threads;
		}
	}

	pthread_mutex_lock(&read_matrix_mutex);
	pthread_cond_wait(&read_matrix_condvar, &read_matrix_mutex);
	pthread_mutex_unlock(&read_matrix_mutex);

    if(inversion_result) {
        fprintf(stderr, "ERROR: matrix is not invertible\n");
        exit_code = 6;
        goto free_threads;
    }

    printf("Inverted matrix:\n");
    print_matrix(inverse, n, n, m);
    printf("\n");

    if(read_matrix(matrix, n, k, filename)) {
        exit_code = 4;
        goto free_args;
    }

	pthread_cond_broadcast(&read_matrix_condvar);

	for(int i = 0; i < threads_amount; i++) {
		pthread_join(threads[i], NULL);
		fprintf(stderr, "kek\n");
		residual_value += args[i].residual_part;
	}

    printf("Residual: %e\n", sqrt(residual_value));
	printf("Total threads time: %.2lf seconds\n",
			(double)thread_total_time / 100);
	printf("Average threads time: %.2lf seconds\n",
			((double)thread_total_time / threads_amount) / 100);

	free_threads:
	free(threads);
	free_args:
	free(args);
    free_inverse:
    free(inverse);
    free_matrix:
    free(matrix);
    final:
    return exit_code;
}

void *thread_execute(void *p_args) {
	long int start_time, finish_time;
	int result;
	struct thread_args *args = (struct thread_args*)p_args;

	start_time = get_thread_time();
	result = invert_matrix(args->matrix, args->inverse_matrix, args->order,
			args->thread_id, args->threads_amount);
	finish_time = get_thread_time();

	pthread_mutex_lock(&total_time_mutex);
	thread_total_time += (finish_time - start_time);
	inversion_result = inversion_result | result;

	if(args->thread_id == 0) {
		pthread_cond_signal(&read_matrix_condvar);
	}
	pthread_mutex_unlock(&total_time_mutex);

	if(result) {
		return NULL;
	}
	
	pthread_mutex_lock(&read_matrix_mutex);
	pthread_cond_wait(&read_matrix_condvar, &read_matrix_mutex);
	pthread_mutex_unlock(&read_matrix_mutex);
	args->residual_part = residual(args->matrix, args->inverse_matrix,
			args->order, args->thread_id, args->threads_amount);

	return NULL;
}
