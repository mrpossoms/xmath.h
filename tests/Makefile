TARGET=$(shell ${CC} -dumpmachine)
INC+=-I../inc
LIBS+=-L../lib/$(TARGET) -L/usr/lib
LINK+=-lm

CFLAGS+=-std=c++11 -g

TESTS=$(patsubst %.c,bin/%,$(wildcard *.c))

bin:
	mkdir bin

bin/%: %.c bin
	$(CXX) $(CFLAGS) $(INC) $(LIBS) $< $(SRC_OBJS) -o $@ $(LINK)

.PHONY: test clean
test: $(TESTS)
	@./.test-runner.py
clean:
	rm -rf bin/
