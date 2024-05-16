CC=mpicc
CFLAGS= -O3 -lm -Isrc -fopenmp -I/usr/lib/x86_64-linux-gnu/openmpi/include/
TEST_FLAGS=-lm -lcunit -fopenmp -Isrc $(shell pkg-config --cflags --libs cunit)

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
	./cg -r benchmark/1_cg_random.in

runlong: cg	
	./cg -r benchmark/cg_large_random.in

# Cleaning up
clean:
	rm -f cg test_cg src/*.o test/*.o
