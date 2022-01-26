CC=mpicc

.PHONY: clean
all: halo

halo: src.c
	$(CC) src.c -lm -o halo

clean:
	rm -f halo
