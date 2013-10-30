EXEC=clpair
#CFLAGS = -Wall -Wextra -I$(INCLUDEDIR) -I$(PYTHIADIR)/include -I$(GSLDIR) -L$(PYTHIADIR)/lib/archive -lpythia8 -llhapdfdummy -O4 -fPIC
CFLAGS = -Wall -Wextra -I$(INCLUDEDIR) -I$(PYTHIADIR)/include -I$(GSLDIR) -L$(PYTHIADIR)/lib/archive -lpythia8 -llhapdfdummy -g
LDFLAGS = -lgsl -lgslcblas $(PYTHIADIR)/lib/archive/libpythia8.a $(PYTHIADIR)/lib/archive/liblhapdfdummy.a
#LDFLAGS = -lgsl -lgslcblas -Linclude/pythia6 -lPythia6
#LDFLAGS = $(gsl-config --libs)
VPATH = src include
############################################
CPP_FILES = $(wildcard src/*.cpp)
HPP_FILES = $(wildcard includes/*.h)
LIB_FILES = $(patsubst src/%.cpp,obj/%.o,$(CPP_FILES))
OBJ_FILES = main.o $(LIB_FILES)
############################################
CC = g++
RM = rm -f
############################################
INCLUDEDIR = include
PYTHIADIR = include/pythia8175
GSLDIR = /usr/include/gsl
############################################

#$(info $(OBJ_FILES))

all: $(EXEC) xsect

$(EXEC): $(OBJ_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: all

obj/%.o: %.cpp %.h
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	$(RM) obj/*.o $(EXEC)

cleanest: clean
	$(RM) mcgen

doc: $(CPP_FILES) $(HPP_FILES) Doxyfile
	doxygen

plot: tryplot.o gnuplot.o
	$(CC) -o $@ $^ $(LDFLAGS)

xsect: utils/xsect.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

probe: utils/probe.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

