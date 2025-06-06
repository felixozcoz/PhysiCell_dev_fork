# versión y nombre del ejecutable
VERSION := $(shell grep . ../../VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := basic_functionality

CC := g++
ifdef PHYSICELL_CPP 
	CC := $(PHYSICELL_CPP)
endif

ifndef STATIC_OPENMP
	STATIC_OPENMP = -fopenmp
endif

ARCH := native

CFLAGS := -march=$(ARCH) -fsanitize=address -O1 -g -mfpmath=both -frtti -fopenmp -m64 -std=c++14

ifeq ($(OS),Windows_NT)
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		UNAME_P := $(shell uname -p)
		var := $(shell which $(CC) | xargs file)
		ifeq ($(lastword $(var)),arm64)
		  CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -fopenmp -m64 -std=c++11
		endif
	endif
endif

CFLAGS_LINK := $(CFLAGS)
COMPILE_COMMAND := $(CC) $(CFLAGS)  $(EXTRA_FLAGS)
LINK_COMMAND := $(CC) $(CFLAGS_LINK) $(EXTRA_FLAGS)

# Objetos de PhysiCell: rutas relativas a la carpeta del makefile (unit_tests/rtree)
BioFVM_OBJECTS := ../../BioFVM/BioFVM_vector.o ../../BioFVM/BioFVM_mesh.o ../../BioFVM/BioFVM_microenvironment.o ../../BioFVM/BioFVM_solvers.o ../../BioFVM/BioFVM_matlab.o \
../../BioFVM/BioFVM_utilities.o ../../BioFVM/BioFVM_basic_agent.o ../../BioFVM/BioFVM_MultiCellDS.o ../../BioFVM/BioFVM_agent_container.o 

PhysiCell_core_OBJECTS := ../../core/PhysiCell_phenotype.o ../../core/PhysiCell_cell_container.o ../../core/PhysiCell_standard_models.o \
../../core/PhysiCell_cell.o ../../core/PhysiCell_custom.o ../../core/PhysiCell_utilities.o ../../core/PhysiCell_constants.o ../../core/PhysiCell_basic_signaling.o \
../../core/PhysiCell_signal_behavior.o ../../core/PhysiCell_rules.o

PhysiCell_module_OBJECTS := ../../modules/PhysiCell_SVG.o ../../modules/PhysiCell_pathology.o ../../modules/PhysiCell_MultiCellDS.o \
../../modules/PhysiCell_various_outputs.o ../../modules/PhysiCell_pugixml.o ../../modules/PhysiCell_settings.o ../../modules/PhysiCell_geometry.o

pugixml_OBJECTS := ../../BioFVM/pugixml.o

ALL_OBJECTS := $(BioFVM_OBJECTS) $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)

# Archivo de test unitario
TEST_SRC := basic_functionality.cpp
TEST_OBJ := basic_functionality.o

.PHONY: all clean

all: $(TEST_OBJ) $(ALL_OBJECTS)
	$(LINK_COMMAND) -o $(PROGRAM_NAME) $(ALL_OBJECTS) $(TEST_OBJ)

$(TEST_OBJ): $(TEST_SRC)
	$(COMPILE_COMMAND) -c $(TEST_SRC) -o $(TEST_OBJ)

# Puedes añadir reglas para compilar los .o del proyecto si quieres, 
# o asumir que ya están compilados en la raíz y sólo enlazar

clean:
	rm -f $(TEST_OBJ) $(PROGRAM_NAME) $(ALL_OBJECTS)
