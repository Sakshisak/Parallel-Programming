CC=mpicc

.PHONY: clean
all: coll

coll: src.c
	$(CC) src.c -lm -o coll

clean:
	rm -f coll
