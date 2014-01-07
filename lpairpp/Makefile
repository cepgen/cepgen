EXEC=clpair
#USE_PYTHIA6 = 1
############################################
CFLAGS = -Wall -Wextra -I$(INCLUDEDIR) -I$(GSLDIR) -llhapdfdummy -g
FFLAGS = -w -g
LDFLAGS = -lgsl -lgslcblas -I$(INCLUDEDIR) -I$(GSLDIR) -lgfortran
VPATH = src include
############################################
CPP_FILES = $(wildcard src/*.cpp)
HPP_FILES = $(wildcard includes/*.h)
LIB_FILES = $(patsubst src/%.cpp,obj/%.o,$(CPP_FILES)) $(patsubst external/%.f,obj/%.o,$(PYTHIA6SRC))
OBJ_FILES = main.o $(LIB_FILES)
############################################
CC = g++
CF = gfortran
RM = rm -f
############################################
INCLUDEDIR = include/
PYTHIA6SRC = $(wildcard external/pythia-6.*.f)
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

obj/%.o: external/%.f
	$(CF) -c $(FFLAGS) $< -o $@

clean:
	$(RM) obj/*.o $(EXEC)

cleanest: clean
	$(RM) mcgen

doc: $(CPP_FILES) $(HPP_FILES) Doxyfile
	doxygen

xsect: utils/xsect.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

probe: utils/probe.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

intest: utils/inelasticparticle.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)
