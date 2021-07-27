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

#include <pthread.h>
#include <sys/resource.h>

#include "common.h"

double f(int n, int k, int i, int j) {
		switch (k) {
		case 1:
			return n - MAX(i, j) + 1;
		case 2:
			return MAX(i, j);
		case 3:
			return ABS(i - j);
		case 4:
			return 1.0/(double)(i + j - 1);
		default:
			return 0;
	}
}

void synchronize(int threads_amount) {
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;
	pthread_mutex_lock(&mutex);
	threads_in++;
	if(threads_in >= threads_amount) {
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} else {
		while(threads_in < threads_amount) {
			pthread_cond_wait(&condvar_in, &mutex);
		}
	}
	threads_out++;
	if(threads_out >= threads_amount) {
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} else {
		while(threads_out < threads_amount) {
			pthread_cond_wait(&condvar_out, &mutex);
		}
	}
	pthread_mutex_unlock(&mutex);
}

long int get_thread_time(void) {
	struct rusage buf;

	getrusage(RUSAGE_SELF, &buf);

	return buf.ru_utime.tv_sec * 100 + buf.ru_utime.tv_usec / 10000;
}
