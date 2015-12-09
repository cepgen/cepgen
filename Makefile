EXEC=clpair
############################################
PYTHIA6SRC = $(wildcard external/pythia-6.*.f)
JETSET7SRC = $(wildcard external/jetset7*.f)
HERWIG6SRC = $(wildcard external/herwig6*.f)
EXTERNALSRC= $(wildcard external/*.f)
INCLUDEDIR = -Iprocesses/ -Iinclude/ -Iexternal/ -Ihadronisers/
############################################
SVNDEV = 'SVN_REV="$(shell svnversion -nq .)"'
CFLAGS     = -Wall -Wextra -fexceptions -Wpointer-arith \
	     $(INCLUDEDIR) -lgsl -g
LDFLAGS    = $(INCLUDEDIR) -lgfortran -lgsl -lgslcblas -Wl,-O2
#LDFLAGS    = $(INCLUDEDIR) -lgfortran -Wl,-O2
FFLAGS     = -w -g
VPATH      = src:include:processes:hadronisers:$(PYTHIA8SRC)/include
############################################
CPP_FILES  = $(wildcard src/*.cpp)
PRO_FILES  = $(wildcard processes/*.cpp)
HAD_FILES  = $(wildcard hadronisers/*.cpp)
HPP_FILES  = $(wildcard include/*.h,external/*.h)
LIB_FILES  = $(patsubst src/%.cpp,obj/%.o,$(CPP_FILES)) \
	     $(patsubst processes/%.cpp,obj/%.o,$(PRO_FILES)) \
	     $(patsubst hadronisers/%.cpp,obj/%.o,$(HAD_FILES)) \
	     $(patsubst external/%.f,obj/%.fo,$(EXTERNALSRC)) 
############################################
CC = @g++
#CC = @clang++
CF = @gfortran
RM = rm -f
############################################

############# FIXME : ROOT #################
RFLAGS = $(shell root-config --cflags)
RLIBS = $(shell root-config --libs)
RHEAD = $(shell root-config --incdir)
############################################

##$(info $(OBJ_FILES))

.PHONY: all

all: $(EXEC)

$(EXEC): main.o $(LIB_FILES)
	@echo "Linking $<..."
	$(CC) -g -o $@ $^ $(LDFLAGS)

diffvm: diffvm.o $(LIB_FILES)
	@echo "Linking $<..."
	$(CC) -g -o $@ $^ $(LDFLAGS)

pptoll: cpptoll.o $(LIB_FILES)
	@echo "Linking $<..."
	$(CC) -g -o $@ $^ $(LDFLAGS)

obj/%.o: %.cpp %.h
	@echo "Building $<..."
	$(CC) -c $(CFLAGS) $< -o $@

obj/%.fo: external/%.f
	@echo "Building (F77) $<..."
	$(CF) -c $(FFLAGS) $< -o $@

obj/%.oxx: %.cxx
	@echo "Building (ROOT) $<..."
	$(CC) -c $(CFLAGS) -I$(RHEAD) $(RFLAGS) $< -o $@

plots/%.oxx: plots/%.cxx
	@echo "Building (ROOT) $<..."
	$(CC) -c $(CFLAGS) -I$(RHEAD) $(RFLAGS) $< -o $@

clean:
	$(RM) obj/*.o $(EXEC) test

doc: $(CPP_FILES) $(HPP_FILES) Doxyfile
	doxygen
	cd doc/latex && make && gnome-open refman.pdf &

xsect: utils/xsect.o $(LIB_FILES)
	@echo "Linking $<..."
	$(CC) -o $@ $^ $(LDFLAGS)

probe: utils/probe.o $(LIB_FILES)
	@echo "Linking $<..."
	$(CC) -o $@ $^ $(LDFLAGS)

intest: utils/inelasticparticle.o $(LIB_FILES)
	@echo "Linking $<..."
	$(CC) -o $@ $^ $(LDFLAGS) -I$(PYTHIA8SRC)/include/

plotter: plots/main.oxx $(LIB_FILES)
	$(CC) -o $@ $^ $(LDFLAGS) $(RLIBS)

test: obj/test.oxx $(LIB_FILES)
	$(CC) -g -o $@ $^ $(LDFLAGS) $(RLIBS)
