EXEC=clpair
#USE_PYTHIA6 = 1
############################################
PYTHIA6SRC = $(wildcard external/pythia-6.*.f)
JETSET7SRC = $(wildcard external/jetset7*.f)
PYTHIA8LIB = $(wildcard external/libpythia8*.so)
INCLUDEDIR = include/
############################################
CFLAGS     = -Wall -Wextra -I$(INCLUDEDIR) -llhapdfdummy -g
FFLAGS     = -w -g
LDFLAGS    = -I$(INCLUDEDIR) -lgfortran
VPATH      = src include
############################################
CPP_FILES  = $(wildcard src/*.cpp)
HPP_FILES  = $(wildcard includes/*.h)
LIB_FILES  = $(patsubst src/%.cpp,obj/%.o,$(CPP_FILES)) \
	     $(patsubst external/%.f,obj/%.fo,$(PYTHIA6SRC)) \
	     $(patsubst external/%.f,obj/%.fo,$(JETSET7SRC))
OBJ_FILES  = main.o $(LIB_FILES)
############################################
CC = g++
CF = gfortran
RM = rm -f
############################################

#$(info $(OBJ_FILES))

all: $(EXEC) xsect

$(EXEC): $(OBJ_FILES)
	$(CC) -g -o $@ $^ $(LDFLAGS)

.PHONY: all

obj/%.o: %.cpp %.h
	$(CC) -c $(CFLAGS) $< -o $@

obj/%.fo: external/%.f
	$(CF) -c $(FFLAGS) $< -o $@

clean:
	$(RM) obj/*.o $(EXEC)

doc: $(CPP_FILES) $(HPP_FILES) Doxyfile
	doxygen
	cd doc/latex && make && gnome-open refman.pdf &

xsect: utils/xsect.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

probe: utils/probe.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

intest: utils/inelasticparticle.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)
