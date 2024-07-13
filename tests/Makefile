TARGET=$(shell ${CC} -dumpmachine)
INC+=-I../inc -I../gitman_sources/asciichart/include/ascii
LIBS+=
LINK+=-lm
CC=gcc
CXXFLAGS+=-std=c++17 -g -O3
CFLAGS=-std=c99 -g -O3

TESTS=$(patsubst %.cpp,bin/%,$(wildcard *.cpp))
TESTS+=$(patsubst %.c,bin/%_c,$(wildcard *.c))

bin:
	mkdir bin

bin/%: %.cpp bin ../inc/xmath.h
	$(CXX) $(CXXFLAGS) $(INC) $(LIBS) $< $(SRC_OBJS) -o $@ $(LINK)

bin/%_c: %.c bin ../inc/xmath.h
	$(CC) $(CFLAGS) $(INC) $(LIBDS) $< -o bin/$<.test $(LINK)

.PHONY: test clean
test: $(C_TESTS) $(TESTS)
	@./.test-runner.py
clean:
	rm -rf bin/
