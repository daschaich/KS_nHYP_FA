CC = gcc
CFLAGS = -std=c99 -O3 -Wall
LIBDIR = ../../gsl/gsl-1.15/lib
INCLUDE_FLAGS = -I../../gsl/gsl-1.15/include
LIBRARY_FLAGS = -lgsl -lgslcblas -lm

main:
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) -c approx_for_modenumber.c
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) -L$(LIBDIR) $(LIBRARY_FLAGS) -o approx_for_modenumber approx_for_modenumber.o

clean:
	rm -f approx_for_modenumber approx_for_modenumber.o
