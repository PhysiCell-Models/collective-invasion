#VERSION := 1.4.1
VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := AMIGOS-invasion

# CC := g++
# CC := g++-mp-7 # typical macports compiler name
CC := g++ # typical homebrew compiler name

# Check for environment definitions of compiler 
# e.g., on CC = g++-7 on OSX
ifdef PHYSICELL_CPP
	CC := $(PHYSICELL_CPP)
endif

ARCH := native # best auto-tuning
# ARCH := core2 # a reasonably safe default for most CPUs since 2007
# ARCH := corei7
# ARCH := corei7-avx # earlier i7 
# ARCH := core-avx-i # i7 ivy bridge or newer 
# ARCH := core-avx2 # i7 with Haswell or newer
# ARCH := nehalem
# ARCH := westmere
# ARCH := sandybridge # circa 2011
# ARCH := ivybridge   # circa 2012
# ARCH := haswell     # circa 2013
# ARCH := broadwell   # circa 2014
# ARCH := skylake     # circa 2015
# ARCH := bonnell
# ARCH := silvermont
# ARCH := skylake-avx512
# ARCH := nocona #64-bit pentium 4 or later 

# CFLAGS := -march=$(ARCH) -Ofast -s -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11
CFLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11
#CFLAGS := -g -march=$(ARCH) -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11

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

COMPILE_COMMAND := $(CC) $(CFLAGS) 

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o 

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o \
PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o PhysiCell_basic_signaling.o \
PhysiCell_signal_behavior.o PhysiCell_rules.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o PhysiCell_geometry.o

# put your custom objects here (they should be in the custom_modules directory)

PhysiCell_custom_module_OBJECTS := AMIGOS-invasion_uncoupled.o extracellular_matrix.o cell_ECM_interactions.o

PhysiCell_custom_module_OBJECT_2 := invasive_spheroid.o extracellular_matrix.o cell_ECM_interactions.o

PhysiCell_custom_module_OBJECT_3 := fibrosis.o extracellular_matrix.o cell_ECM_interactions.o

pugixml_OBJECTS := pugixml.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS)  $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)

LEADER_FOLLOWER_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS)
INVASIVE_SPHEROID_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECT_2)
FIBROSIS_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECT_3)
#compile the project 
	
all: main-ecm.cpp $(LEADER_FOLLOWER_OBJECTS)
	$(COMPILE_COMMAND) -o $(PROGRAM_NAME) $(LEADER_FOLLOWER_OBJECTS) main-ecm.cpp 

invasive_spheroid: main_invasive_spheroid.cpp $(INVASIVE_SPHEROID_OBJECTS)
	$(COMPILE_COMMAND) -o invasive_spheroid $(INVASIVE_SPHEROID_OBJECTS) main_invasive_spheroid.cpp 

fibrosis: main_fibrosis.cpp $(FIBROSIS_OBJECTS)
	$(COMPILE_COMMAND) -o fibrosis $(FIBROSIS_OBJECTS) main_fibrosis.cpp 

# PhysiCell core components	

PhysiCell_phenotype.o: ./core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_phenotype.cpp
	
PhysiCell_digital_cell_line.o: ./core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: ./core/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell.cpp 

PhysiCell_cell_container.o: ./core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell_container.cpp 
	
PhysiCell_standard_models.o: ./core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_standard_models.cpp 
	
PhysiCell_utilities.o: ./core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_utilities.cpp 
	
PhysiCell_custom.o: ./core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_custom.cpp 
	
PhysiCell_constants.o: ./core/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_constants.cpp 

PhysiCell_signal_behavior.o: ./core/PhysiCell_signal_behavior.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_signal_behavior.cpp 

PhysiCell_rules.o: ./core/PhysiCell_rules.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_rules.cpp 

# BioFVM core components (needed by PhysiCell)
	
BioFVM_vector.o: ./BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_vector.cpp 

BioFVM_agent_container.o: ./BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_agent_container.cpp 
	
BioFVM_mesh.o: ./BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_mesh.cpp 

BioFVM_microenvironment.o: ./BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_microenvironment.cpp 

BioFVM_solvers.o: ./BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_solvers.cpp 

BioFVM_utilities.o: ./BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_utilities.cpp 
	
BioFVM_basic_agent.o: ./BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_basic_agent.cpp 
	
BioFVM_matlab.o: ./BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: ./BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_MultiCellDS.cpp
	
pugixml.o: ./BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/pugixml.cpp
	
# standard PhysiCell modules

PhysiCell_SVG.o: ./modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: ./modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: ./modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: ./modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_various_outputs.cpp
	
PhysiCell_pugixml.o: ./modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pugixml.cpp
	
PhysiCell_settings.o: ./modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_settings.cpp	

PhysiCell_basic_signaling.o: ./core/PhysiCell_basic_signaling.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_basic_signaling.cpp 	
	
PhysiCell_geometry.o: ./modules/PhysiCell_geometry.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_geometry.cpp 
	
# user-defined PhysiCell modules

extracellular_matrix.o: ./custom_modules/extracellular_matrix.cpp 
	$(COMPILE_COMMAND) -c ./custom_modules/extracellular_matrix.cpp

cell_ECM_interactions.o: ./custom_modules/cell_ECM_interactions.cpp 
	$(COMPILE_COMMAND) -c ./custom_modules/cell_ECM_interactions.cpp

AMIGOS-invasion_uncoupled.o: ./custom_modules/AMIGOS-invasion_uncoupled.cpp 
	$(COMPILE_COMMAND) -c ./custom_modules/AMIGOS-invasion_uncoupled.cpp

invasive_spheroid.o: ./custom_modules/invasive_spheroid.cpp
	$(COMPILE_COMMAND) -c ./custom_modules/invasive_spheroid.cpp

fibrosis.o: ./custom_modules/fibrosis.cpp
	$(COMPILE_COMMAND) -c ./custom_modules/fibrosis.cpp

# cleanup

#reset:
#	rm -f *.cpp 
#	cp ./sample_projects/Makefile-default Makefile 
#	rm -f ./custom_modules/*
#	touch ./custom_modules/empty.txt 
	
clean:
	rm -f *.o
	rm -f $(PROGRAM_NAME)*

data-cleanup:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
	rm -f ./output/*
	touch ./output/empty.txt

data-cleanup-all:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
	rm -f *.pov
	rm -f ./Output/*
	rm -f ./SVG/*

data-cleanup-light:
	rm -f *.mat
	rm -f *.xml
	rm -f *.svg
	rm -f *.pov
	rm -f ./Output/*.mat
	rm -f ./Output/*.xml
	rm -f ./Output/*.svg
	rm -f ./Output/*.png
	rm -f ./SVG/*
	touch ./output/empty.txt

# archival 
	
zip:
	zip -r latest.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/* 
	cp latest.zip $$(date +%b_%d_%Y_%H%M).zip
	cp latest.zip VERSION_$(VERSION).zip 
	mv *.zip archives/

zip-data:
	zip -r latest-data.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* output/*.xml output/*.mat output/*.svg output/*.png output/*.m output/*.mov output/*.gif output/*.mp4
	cp latest-data.zip $$(date +%b_%d_%Y_%H%M).zip
	cp latest-data.zip VERSION_$(VERSION).zip 
	mv *.zip data_archives/

zip-source:
	zip -r latest-data.zip Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* python_imaging/* modules/* output/*.m
	cp latest-data.zip $$(date +%b_%d_%Y_%H%M).zip
	cp latest-data.zip VERSION_$(VERSION).zip 
	mv *.zip data_archives/

move-data:
	mkdir $$(date +%b_%d_%Y_%H%M)
	cp -a Makefile* /$$(date +%b_%d_%Y_%H%M)
	mv /$$(date +%b_%d_%Y_%H%M) data//$$(date +%b_%d_%Y_%H%M)
	
tar:
	tar --ignore-failed-read -czf latest.tar Makefile* *.cpp *.h BioFVM/* config/* core/* custom_modules/* matlab/* modules/* sample_projects/* 
	cp latest.tar $$(date +%b_%d_%Y_%H%M).tar
	cp latest.tar VERSION_$(VERSION).tar
	mv *.tar archives/

unzip: 
	cp ./archives/latest.zip . 
	unzip latest.zip 
	
untar: 
	cp ./archives/latest.tar .
	tar -xzf latest.tar
