#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <sched.h>
#include <omp.h>

double gettime(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ((double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec);
}

int main(int argc, char **argv)
{
    if (argc < 2)
        return 1;
    double start, total;
	double pi, sum = 0.0;
    long long num_steps = 1000 * 1000 * 1000;
	double step = 1. / (double)num_steps;
    int threads_num;
    double *pi_tab;
    int n;
    int i;

    int opt = atoi(argv[1]);
    switch (opt)
    {
        case 1:
            /* omp_set_num_threads(4); */
            start = gettime();
#pragma omp parallel for reduction (+:sum)
            for (int i = 0; i < num_steps; i++)
            {
                double x = (i + 0.5)*step;
                sum += 4.0 / (1. + x*x);
            }
            pi = sum * step;
            total = gettime() - start;
            break;

        case 2:
            threads_num = 4;
            omp_set_num_threads(threads_num);
            pi_tab = malloc(threads_num * sizeof(double));
            for (int i = 0; i < threads_num; i++)
                pi_tab[i] = 0;
            start = gettime();
#pragma omp parallel
            {
                int ind = omp_get_thread_num();
#pragma omp for
                for (int i = 0; i < num_steps; i++)
                {

                    double x = (i + .5) * step;
                    pi_tab[ind] += 4.0 / (1. + x*x);
                }
            }
            for (int i = 0; i < threads_num; i++)
                sum += pi_tab[i];
            pi = sum * step;
            total = gettime() - start;
            free(pi_tab);
            break;

        case 3:
            threads_num = 2;
            omp_set_num_threads(threads_num);
            n = 20;
            pi_tab = malloc(n * sizeof(double));
            for (int j = 0; j < n-1; j++)
            {
                pi_tab[j] = 0;
                pi_tab[j + 1] = 0;
                start = gettime();
#pragma omp parallel 
                {
                    int ind = omp_get_thread_num();
                    /* SetThreadAffinityMask(GetCurrentThread(), 1 << omp_get_thread_num()+1); */
                    cpu_set_t mask;
                    CPU_ZERO(&mask);
                    CPU_SET(ind & 1, &mask);
                    sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set_t), &mask);
#pragma omp for
                    for (int i = 0; i < num_steps; i++)
                    {
                        double x = (i + .5)*step;
                        pi_tab[j+ind] += 4.0 / (1. + x*x);
                    }
                }
                sum = pi_tab[j] + pi_tab[j+1];		
                pi = sum*step;
                total = gettime() - start;
                printf("%d:", j);
                printf("Wartosc liczby PI wynosi %15.12f\n", pi);
                printf("Czas przetwarzania wynosi %lf sekund\n", total);
            }
            free(pi_tab);
            break;

        case 4:
            omp_set_num_threads(4);
            start = gettime();
#pragma omp parallel  reduction (+:sum)
            {
                /* SetThreadAffinityMask(GetCurrentThread(), 1 << omp_get_thread_num()); */
#pragma omp for
                for (int i = 0; i < num_steps; i++)
                {
                    double x = (i + 0.5)*step;
                    sum += 4.0 / (1. + x*x);
                }
            }
            pi = sum*step;
            total = gettime() - start;
            break;

        case 5:
            start = gettime();
#pragma omp parallel  reduction (+:sum)
            {
                /* SetThreadAffinityMask(GetCurrentThread(), 1 << (omp_get_thread_num()%2 +1)); */
#pragma omp for
                for (int i = 0; i < num_steps; i++)
                {
                    double x = (i + 0.5)*step;
                    sum += 4.0 / (1. + x*x);
                }
            }
            pi = sum * step;
            total = gettime() - start;
            break;
    }
    printf("Wartosc liczby PI wynosi %15.12lf\n", pi);
	printf("Czas przetwarzania wynosi %lf sekund\n", total);
}
