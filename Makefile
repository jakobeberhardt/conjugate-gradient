CC=gcc
CFLAGS=-lm -Isrc
TEST_FLAGS=-lm -lcunit -Isrc $(shell pkg-config --cflags --libs cunit)

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

# Cleaning up
clean:
	rm -f cg test_cg src/*.o test/*.o
