
CFLAGS=-Wall -Wextra --pedantic -O3
# CFLAGS=-Wall -Wextra --pedantic -g -fsanitize=undefined
# CFLAGS=-Wall -Wextra --pedantic -g -fsanitize=address

all:	OpenSimplex2F.o test_OpenSimplex2F

OpenSimplex2F.o:	OpenSimplex2F.c OpenSimplex2F.h Makefile
	$(CC) ${CFLAGS} -c OpenSimplex2F.c -lm

test_OpenSimplex2F:	test_OpenSimplex2F.c OpenSimplex2F.o Makefile
	$(CC) ${CFLAGS} -o test_OpenSimplex2F test_OpenSimplex2F.c OpenSimplex2F.o -lm -lpng

clean:
	rm -f *.o test_OpenSimplex2F test2d.png test3d.png test4d.png


