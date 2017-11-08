ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXXFLAGS      = -g -Wall $(ROOTCFLAGS)

CC     = gcc
CXX    = c++ -w
SRC    = ./src
EXE    = ./exe
BIN    = ./bin
BUILD  = ./build
LIB    = ./lib
INC    = ./inc
ADL    = ../ADL4

LIBS    = -lm  $(ROOTLIBS) -lTreePlayer -lRGL  $(SYSLIBS)
LIBS   += -L${ADL}/lib/ -lADL-4-2
PROG    = SimulatePulse GetPulserResponse ConvolutePulses

INCFLAGS      = -I$(INC)/ -I${ADL}/include/ -I/lfs/l3/gerda/kermaidy/Analysis/software/src/gelatio/Decoders/
CXXFLAGS     += $(INCFLAGS)
LOCALLIBS     = -L$(LIB) -lSIMPULSE

# MGDO
CXXFLAGS += -I/lfs/l3/gerda/kermaidy/Analysis/software/gerda/linux-scientific-7.2-x86_64/master/include/mgdo
LIBS += -L/lfs/l3/gerda/kermaidy/Analysis/software/gerda/linux-scientific-7.2-x86_64/master/lib -lMGDOBase -lMGDOTransforms -lMGDORoot -lMGDOGerda

EXESRC = $(wildcard $(EXE)/*.cc)
EXEOBJ = $(subst $(EXE),$(BIN)/,$(patsubst %.cc,%_x,$(EXESRC)))

SRCSRC = $(wildcard $(SRC)/*.cc)
SRCOBJ = $(subst $(SRC)/,$(BUILD)/,$(patsubst %.cc,%.o,$(SRCSRC)))

default: all
all: tools library exe
	@echo "linking executable"
	@ln -f -s bin/SimulatePulse_x ./SimulatePulse
	@ln -f -s bin/GetPulserResponse_x ./GetPulserResponse
	@ln -f -s bin/ConvolutePulses_x ./ConvolutePulses

tools : $(SRCOBJ)
exe:    $(EXEOBJ)

library:
	@echo "Building library libSIMPULSE.so"
	@$(CXX) -shared -o $(LIB)/libSIMPULSE.so $(BUILD)/*.o

$(BUILD)/%.o: $(SRC)/%.cc $(INC)/%.hh
	@echo "in ["$(SRC)"] -> Building object : " $@
	@$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $<

$(BIN)/%_x: $(EXE)/%.cc
	@echo "Building executable : " $@
	@$(CXX) $(CXXFLAGS) -o $@ $?  $(LOCALLIBS) $(LIBS)

clean:
	rm -f SimulatePulse bin/* build/* lib/* JobErrors/* JobOutput/*

distclean: clean

