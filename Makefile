test: test.c
	gcc -s -std=c99 -fopenmp -D_GNU_SOURCE -D_BSD_SOURCE -D_POSIX_C_SOURCE=200809L test.c -o test
