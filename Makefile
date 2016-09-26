EXEC=cepgen
############################################
PYTHIA6SRC = $(wildcard external/pythia-6.*.f)
JETSET7SRC = $(wildcard external/jetset7*.f)
HERWIG6SRC = $(wildcard external/herwig6*.f)
EXTERNALSRC= $(wildcard external/*.f)
INCLUDEDIR = -Iprocesses/ -Iinclude/ -Iexternal/ -Ihadronisers/ -Iexport/
OBJDIR     = obj/
VPATH      = src:include:processes:hadronisers:export
ifdef PYTHIA8
  $(info Using Pythia 8 as included from PYTHIA8=$(PYTHIA8))
  #make -f $(PYTHIA8)/Makefile
  INCLUDEDIR += -I$(PYTHIA8)/include/ -L$(PYTHIA8)/lib/archive/ -lpythia8 -llhapdfdummy -DPYTHIA8=1
else
  #$(info PYTHIA8 variable is not set... skipping its compilation)
endif
############################################
CFLAGS     = -Wall -Wextra -fexceptions -Wpointer-arith \
	     $(INCLUDEDIR) -g -pedantic-errors
LDFLAGS    = $(INCLUDEDIR) -lgfortran -lgsl -lgslcblas -Wl,-O2
FFLAGS     = -w -g
HEPMCFLAGS = -I/usr/local/include/HepMC/ -L/usr/local/lib64/ -lHepMC
#HEPMCFLAGS = -lHepMC
#INCLUDEDIR += -I/usr/local/include/HepMC/
############################################
CPP_FILES  = $(wildcard src/*.cpp)
PRO_FILES  = $(wildcard processes/*.cpp)
HAD_FILES  = $(wildcard hadronisers/Pythia6Hadroniser.cpp)
HPP_FILES  = $(wildcard include/*.h,external/*.h)
LIB_FILES  = $(patsubst src/%.cpp,$(OBJDIR)%.o,$(CPP_FILES)) \
	     $(patsubst processes/%.cpp,$(OBJDIR)%.o,$(PRO_FILES)) \
	     $(patsubst hadronisers/%.cpp,$(OBJDIR)%.o,$(HAD_FILES)) \
	     $(patsubst external/%.f,$(OBJDIR)%.fo,$(EXTERNALSRC))
EXP_FILES  = $(wildcard export/*.cpp)
EXP_LIB_FILES = $(patsubst export/%.cpp,$(OBJDIR)%.o,$(EXP_FILES))
############################################
CC = @g++
#CC = @clang++
CF = @gfortran
RM = rm -f
############################################

############# FIXME : ROOT #################
RFLAGS = $(shell root-config --cflags) -std=c++11
RLIBS = $(shell root-config --libs)
RHEAD = $(shell root-config --incdir)
############################################

##$(info $(OBJ_FILES))

.PHONY: all

all: $(EXEC)

$(EXEC): $(OBJDIR)main.opp $(LIB_FILES)
	@echo "Linking $<..."
	@$(CC) -g -o $@ $^ $(LDFLAGS)

cepgen-lhe: $(OBJDIR)cepgen-lhe.opp $(LIB_FILES) $(EXP_LIB_FILES)
	@echo "Linking $<..."
	@$(CC) -g -o $@ $^ $(LDFLAGS) $(HEPMCFLAGS)

$(OBJDIR)%.o: %.cpp %.h | $(OBJDIR)
	@echo "Building $<..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)%.fo: external/%.f | $(OBJDIR)
	@echo "Building (F77) $<..."
	@$(CF) -c $(FFLAGS) $< -o $@
$(OBJDIR)%.opp: %.cpp | $(OBJDIR)
	#@echo $(LIB_FILES)
	@echo "Building $<..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)%.oxx: %.cxx | $(OBJDIR)
	@echo "Building (ROOT) $<..."
	@$(CC) -c $(CFLAGS) -I$(RHEAD) $(RFLAGS) $< -o $@

plots/%.oxx: plots/%.cxx
	@echo "Building (ROOT) $<..."
	@$(CC) -c $(CFLAGS) -I$(RHEAD) $(RFLAGS) $< -o $@

clean:
	@$(RM) -r $(OBJDIR)
	@$(RM) $(EXEC) test

doc: $(CPP_FILES) $(HPP_FILES) Doxyfile
	doxygen
	cd doc/latex && make && gnome-open refman.pdf &

xsect: utils/xsect.o $(LIB_FILES)
	@echo "Linking $<..."
	@$(CC) -o $@ $^ $(LDFLAGS)

probe: utils/probe.o $(LIB_FILES)
	@echo "Linking $<..."
	@$(CC) -o $@ $^ $(LDFLAGS)

intest: utils/inelasticparticle.o $(LIB_FILES) $(EXP_LIB_FILES)
	@echo "Linking $<..."
	#@$(CC) -o $@ $^ $(LDFLAGS) -I$(PYTHIA8SRC)/include/
	@$(CC) -o $@ $^ $(LDFLAGS) $(HEPMCFLAGS)

plotter: plots/main.oxx $(LIB_FILES)
	@$(CC) -o $@ $^ $(LDFLAGS) $(RLIBS)

test: $(OBJDIR)/test.oxx $(LIB_FILES)
	@$(CC) -g -o $@ $^ $(LDFLAGS) $(RLIBS)

$(OBJDIR):
	mkdir -p $(OBJDIR)
