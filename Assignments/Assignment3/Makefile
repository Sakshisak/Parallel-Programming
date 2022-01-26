CC=mpicc

.PHONY: clean
all: a.out

a.out: src.c
	$(CC) $^ -lm 

clean:
	rm -f a.out
