CC=icc
CFLAGS=-Wall -O3 -vec_report3 -g -lm -openmp

#test_prog21r : HW3ParallelCorrelationFunc1_0.c

SDSS_CF_Parallel : HW3ParallelCorrelationFunc2_2.c
	$(CC) -o $@ $< $(CFLAGS)

.PHONY: clean

clean :
	rm -f SDSS_CF_Parallel



