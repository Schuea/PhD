# Where to find user code.
USER_DIR = ../samples

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include -g -Wall -Wextra -pthread -O2 -std=c++1y $(shell root-config --cflags) -Wfatal-errors 

# Flags passed to the C++ compiler.
CXXFLAGS +=  $(shell root-config --glibs)

# Executables produced by this Makefile
RUN = ../NewOccupancyPlotter ../CompareOccupancy ../EnergyPlotter ../CompareEnergy ../HitTimePlotter ../CompareHitTime ../HitsPerSubdetectorPlotter
#RUN = ../NewOccupancyPlotter ../NewOccupancyPlotter_Phi ../CompareOccupancies

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = #Test

# Useful variables
MAIN_OCCUPANCY_NEW = ../objects/main_occupancy_new.o
MAIN_COMPAREOCCUPANCY = ../objects/main_CompareOccupancy.o
MAIN_HITTIME = ../objects/main_hittime_distribution.o
MAIN_COMPAREHITTIME = ../objects/main_CompareHitTime.o
MAIN_ENERGY = ../objects/main_energy_distribution.o
MAIN_COMPAREENERGY = ../objects/main_CompareEnergy.o
MAIN_HITS = ../objects/Hits_per_Subdetector_MuonSpoilerStudy.o
#MAIN_OCCUPANCY_NEW_PHI = ../objects/main_occupancy_r_phi_new.o
#MAIN_OCCUPANCY_NEW_COMPARISON = ../objects/main_CompareOccupancy.o
USEFUL = ../objects/UsefulFunctions.o
STYLE = ../objects/Style.o
STYLEUTILS = ../objects/Style_Utils.o
GENERALFUNCTIONS_NEW = ../objects/GeneralFunctions_SiDBkgSim_new.o
SUBDETECTOR_NEW = ../objects/Subdetector_class_new.o
CELLHITS_NEW = ../objects/CellHits_class_new.o

INCLUDES = -I ../include/ -I ../test/include/ -I ../../General/include/

OBJECTS = $(CONFIG) $(USEFUL) $(STYLE) $(STYLEUTILS) $(GENERALFUNCTIONS_NEW) $(SUBDETECTOR_NEW) $(CELLHITS_NEW)

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# House-keeping build targets.

all : $(RUN) $(TESTS)

clean :
	@rm -f ./$(RUN) ../objects/*

# Builds gtest.a and gtest_main.a.
# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
../objects/gtest-all.o : $(GTEST_SRCS_)
	@$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@

../objects/gtest_main.o : $(GTEST_SRCS_)
	@$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc -o $@

../objects/gtest.a : ../objects/gtest-all.o
	@$(AR) $(ARFLAGS) $@ $^

../objects/gtest_main.a : ../objects/gtest-all.o ../objects/gtest_main.o
	@$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.

../objects/test.o : ../test/src/test.cc $(GTEST_HEADERS)
	@$(CXX) $(CPPFLAGS) -c ../test/src/test.cc $(INCLUDES) -o $@ $(CXXFLAGS)

Test : ../objects/test.o ../objects/gtest.a $(OBJECTS)
	@$(CXX) $(CPPFLAGS) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

# Building the graviton code

../NewOccupancyPlotter : $(MAIN_OCCUPANCY_NEW) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../CompareOccupancy : $(MAIN_COMPAREOCCUPANCY) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../HitTimePlotter : $(MAIN_HITTIME) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../CompareHitTime : $(MAIN_COMPAREHITTIME) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../EnergyPlotter : $(MAIN_ENERGY) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../CompareEnergy : $(MAIN_COMPAREENERGY) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../HitsPerSubdetectorPlotter : $(MAIN_HITS) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

#../NewOccupancyPlotter_Phi : $(MAIN_OCCUPANCY_NEW_PHI) $(OBJECTS)
#	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

#../CompareOccupancies : $(MAIN_OCCUPANCY_NEW_COMPARISON) $(OBJECTS)
#	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../objects/main_occupancy_new.o : ../src/main_occupancy_new.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_CompareOccupancy.o : ../src/main_CompareOccupancy.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_hittime_distribution.o : ../src/main_hittime_distribution.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_CompareHitTime.o : ../src/main_CompareHitTime.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_energy_distribution.o : ../src/main_energy_distribution.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_CompareEnergy.o : ../src/main_CompareEnergy.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Hits_per_Subdetector_MuonSpoilerStudy.o : ../src/Hits_per_Subdetector_MuonSpoilerStudy.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

#../objects/main_occupancy_r_phi_new.o : ../src/main_occupancy_r_phi_new.cc
#	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

#../objects/main_CompareOccupancy.o : ../src/main_CompareOccupancy.cc
#	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/UsefulFunctions.o : ../../General/src/UsefulFunctions.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Style.o : ../../General/src/Style.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Style_Utils.o : ../../General/src/Style_Utils.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/GeneralFunctions_SiDBkgSim_new.o : ../src/GeneralFunctions_SiDBkgSim_new.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Subdetector_class_new.o : ../src/Subdetector_class_new.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/CellHits_class_new.o : ../src/CellHits_class_new.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@
