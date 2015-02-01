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

pid_t GetCurrentThread()
{
    return syscall(SYS_gettid);
}

void SetThreadAffinityMask(pid_t tid, int cpu)
{
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(cpu, &mask);
    sched_setaffinity(tid, sizeof(cpu_set_t), &mask);
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stderr, "usage: %s CASE NUM_STEPS\n", argv[0]);
        return 1;
    }
    double start, total;
	double pi, sum = 0.0;
    long long num_steps = atoll(argv[2]);
	double step = 1. / (double)num_steps;
    int threads_num;
    double *pi_tab;
    int n;
    int i;

    int opt = atoi(argv[1]);
    switch (opt)
    {
        case 0:
            start = gettime();
            for (int i = 0; i < num_steps; i++)
            {
                double x = (i + 0.5)*step;
                sum += 4.0 / (1. + x*x);
            }
            pi = sum * step;
            total = gettime() - start;
            break;
            
        case 1:
            /* omp_set_num_threads(4); */
            start = gettime();
#pragma omp parallel for
            for (int i = 0; i < num_steps; i++)
            {
                double x = (i + 0.5)*step;
                sum += 4.0 / (1. + x*x);
            }
            pi = sum * step;
            total = gettime() - start;
            break;

        case 2:
            start = gettime();
#pragma omp parallel for
            for (int i = 0; i < num_steps; i++)
            {
                double x = (i + .5) * step;
#pragma omp atomic
                sum += 4.0 / (1. + x*x);
            }
            pi = sum * step;
            total = gettime() - start;
            break;

        case 3:
            start = gettime();
#pragma omp parallel  reduction (+:sum)
            {
                /* SetThreadAffinityMask(GetCurrentThread(), omp_get_thread_num() & 1); */
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

        case 4:
            threads_num = 2;
            omp_set_num_threads(threads_num);
            n = 20;
            pi_tab = (double *)malloc(n * sizeof(double));
            for (int j = 0; j < n-1; j++)
            {
                pi_tab[j] = 0;
                pi_tab[j + 1] = 0;
                start = gettime();
#pragma omp parallel 
                {
                    int ind = omp_get_thread_num();
                    /* SetThreadAffinityMask(GetCurrentThread(), ind & 1); */
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

        case 5:
            start = gettime();
#pragma omp parallel  reduction (+:sum)
            {
                /* SetThreadAffinityMask(GetCurrentThread(), 1); */
                SetThreadAffinityMask(GetCurrentThread(), (omp_get_thread_num() & 1));
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
