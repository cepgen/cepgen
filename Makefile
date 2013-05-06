CC = g++
CFLAGS = -I include/ -Wall -Wextra -I /usr/include/gsl -O4 -fPIC -fstack-usage #-std=c++11 #-dD
LDFLAGS = -lgsl -lgslcblas
VPATH = src include
CPP_FILES = $(wildcard src/*.cpp)
HPP_FILES = $(wildcard includes/*.h)
OBJ_FILES = main.o $(patsubst src/%.cpp,%.o,$(CPP_FILES))
RM = rm -f

$(info $(OBJ_FILES))

all: mcgen clean

mcgen: $(OBJ_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

#mcgen.o: main.cpp mcgen.cpp mcgen.h
#	$(CC) -c $(CFLAGS) $<

.PHONY: all

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	$(RM) *.o

cleanest: clean
	$(RM) mcgen

doc: $(CPP_FILES) $(HPP_FILES) Doxyfile
	doxygen

test: tryplot.o gnuplot.o
	$(CC) -o $@ $^ $(LDFLAGS)

