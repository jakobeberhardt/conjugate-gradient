CC=mpicc
CFLAGS= -O3 -lm -mavx2 -Isrc -fopenmp -I/usr/lib/x86_64-linux-gnu/openmpi/include/
TEST_FLAGS=-lm -lcunit -mavx2 -fopenmp -Isrc $(shell pkg-config --cflags --libs cunit)

all: cg

cg: src/cg.o src/main.o src/misc.o
	$(CC) -o $@ $^ $(CFLAGS)

src/%.o: src/%.c
	$(CC) -c $< -o $@ $(CFLAGS)

test: cg test_cg
	./test_cg

test_cg: test/test_cg.o src/cg.o src/misc.o
	$(CC) -o $@ $^ $(TEST_FLAGS)

test/test_cg.o: test/test_cg.c src/cg.h src/misc.h
	$(CC) -c $< -o $@ $(TEST_FLAGS)

run: cg
	OMP_NUM_THREADS=1 mpirun -np 1 ./cg -r benchmark/3_cg_random.in

runlong: cg	
	./cg -r benchmark/cg_large_random.in

valgrind: cg	
	OMP_NUM_THREADS=1 mpirun -np 1 valgrind --tool=callgrind ./cg -r benchmark/3_cg_random.in && rm callgrind.out.*

clean:
	rm -f cg test_cg src/*.o test/*.o callgrind.out.*
