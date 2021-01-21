TARGET=$(shell ${CC} -dumpmachine)
INC+=-I../inc
LIBS+=
LINK+=-lm
CC=gcc
CXXFLAGS+=-std=c++11 -g
CFLAGS=-std=c99 -g

TESTS=$(patsubst %.cpp,bin/%,$(wildcard *.cpp))
TESTS+=$(patsubst %.c,bin/%_c,$(wildcard *.c))

bin:
	mkdir bin

bin/%: %.cpp bin
	$(CXX) $(CPPFLAGS) $(INC) $(LIBS) $< $(SRC_OBJS) -o $@ $(LINK)

bin/%_c: %.c bin
	$(CC) $(CFLAGS) $(INC) $(LIBDS) $< -o bin/$<.test $(LINK) 

.PHONY: test clean 
test: $(C_TESTS) $(TESTS)
	@./.test-runner.py
clean:
	rm -rf bin/
