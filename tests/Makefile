TARGET=$(shell ${CC} -dumpmachine)
INC+=-I../inc
LIBS+=-L../lib/$(TARGET) -L/usr/lib
LINK+=-lm

CFLAGS+=-std=c++11 -g

TESTS=$(patsubst %.cpp,bin/%,$(wildcard *.cpp))

bin:
	mkdir bin

bin/%: %.cpp bin
	$(CXX) $(CFLAGS) $(INC) $(LIBS) $< $(SRC_OBJS) -o $@ $(LINK)

.PHONY: test clean
test: $(TESTS)
	@./.test-runner.py
clean:
	rm -rf bin/
