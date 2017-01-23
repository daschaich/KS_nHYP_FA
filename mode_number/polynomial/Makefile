CC = gcc
CFLAGS = -std=c99 -O3 -Wall
LIBRARY_FLAGS = -lgsl -lgslcblas -lm

main:
	$(CC) $(CFLAGS) -c approx_for_modenumber.c
	$(CC) $(CFLAGS) $(LIBRARY_FLAGS) -o approx_for_modenumber approx_for_modenumber.o

clean:
	rm -f approx_for_modenumber approx_for_modenumber.o
