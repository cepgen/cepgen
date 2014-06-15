EXEC=clpair
############################################
PYTHIA6SRC = $(wildcard external/pythia-6.*.f)
JETSET7SRC = $(wildcard external/jetset7*.f)
PYTHIA8SRC = $(wildcard external/pythia8*/)
HERWIG6SRC = $(wildcard external/herwig6*.f)
PYTHIA8LIB = $(PYTHIA8SRC)lib
INCLUDEDIR = -Imodels/ -Iinclude/
############################################
CFLAGS     = -Wall -Wextra -fno-exceptions -Wpointer-arith \
	     $(INCLUDEDIR) -D$(SVNDEV) -g #-I$(PYTHIA8SRC)/include
FFLAGS     = -w -g
LDFLAGS    = $(INCLUDEDIR) -lgfortran -L$(PYTHIA8LIB) -I$(PYTHIA8SRC)include -lpythia8 -llhapdfdummy -Wl,-O2
VPATH      = src:include:models:$(PYTHIA8SRC)/include
SVNDEV = 'SVN_REV="$(shell svnversion -nq .)"'
############################################
CPP_FILES  = $(wildcard src/*.cpp)
MOD_FILES  = $(wildcard models/*.cpp)
HPP_FILES  = $(wildcard include/*.h)
LIB_FILES  = $(patsubst src/%.cpp,obj/%.o,$(CPP_FILES)) \
	     $(patsubst models/%.cpp,obj/%.o,$(MOD_FILES)) \
	     $(patsubst external/%.f,obj/%.fo,$(PYTHIA6SRC)) \
	     $(patsubst external/%.f,obj/%.fo,$(HERWIG6SRC)) \
	     $(patsubst external/%.f,obj/%.fo,$(JETSET7SRC)) 
	     #-L$(PYTHIA8LIB) -I$(PYTHIA8SRC)/include -lpythia8
############################################
CC = g++
#CC = clang++
CF = gfortran
RM = rm -f
############################################

############# FIXME : ROOT #################
RFLAGS = $(shell root-config --cflags)
RLIBS = $(shell root-config --libs)
RHEAD = $(shell root-config --incdir)
############################################

##$(info $(OBJ_FILES))

.PHONY: all

all: $(EXEC) #xsect

$(EXEC): main.o $(LIB_FILES)
	$(CC) -g -o $@ $^ $(LDFLAGS) -L$(PYTHIA8LIB) -lpythia8 -llhapdfdummy

#obj/pythia8hadroniser.o: pythia8hadroniser.cpp pythia8hadroniser.h
#	$(CC) -c $(CFLAGS) -llhapdfdummy $< -o $@

obj/%.o: %.cpp %.h
	$(CC) -c $(CFLAGS) $< -o $@

obj/%.fo: external/%.f
	$(CF) -c $(FFLAGS) $< -o $@

obj/%.oxx: %.cxx
	$(CC) -c $(CFLAGS) -I$(RHEAD) $(RFLAGS) $< -o $@

plots/%.oxx: plots/%.cxx
	$(CC) -c $(CFLAGS) -I$(RHEAD) $(RFLAGS) $< -o $@

clean:
	$(RM) obj/*.o $(EXEC) test

doc: $(CPP_FILES) $(HPP_FILES) Doxyfile
	doxygen
	cd doc/latex && make && gnome-open refman.pdf &

xsect: utils/xsect.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

probe: utils/probe.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS)

intest: utils/inelasticparticle.o $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS) -I$(PYTHIA8SRC)/include/

plotter: plots/main.oxx $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS) $(RLIBS)

test: obj/test.oxx $(LIB_FILES)
	$(CC) -g -o $@ $^ $(LDFLAGS) $(RLIBS)
