# Where to find user code.
USER_DIR = ../samples

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include -g -Wall -Wextra -pthread -O2 -std=c++1y $(shell root-config --cflags) -Wfatal-errors 

# Flags passed to the C++ compiler.
CXXFLAGS +=  $(shell root-config --glibs)

# Executables produced by this Makefile
RUN = ../BDSIM_ParticleCount_Plotter ../BDSIM_xy_PrimariesSecondaries_Plotter ../RHULcherenkov_AverageSignal_Aperture_Plotter ../RHULcherenkov_AverageSignal_AsymmAperture_Plotter ../RHULcherenkov_AverageSignal_HalfAperture_Plotter ../RHULcherenkov_AverageSignal_Aperture_VacuumPressures_Plotter ../RHULcherenkov_AverageSignal_VoltageNormalization_Plotter ../RHULcherenkov_AverageNoise_Voltage_Plotter


# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = #Test

# Useful variables
MAIN_OCCUPANCY = ../objects/main_occupancy.o
MAIN_PRIMESECOND = ../objects/main_BDSIM_ParticleCount.o
MAIN_XY_PRIMESECOND = ../objects/main_BDSIM_xy_PrimariesSecondaries.o
MAIN_RHUL_AVERSIGNAL = ../objects/main_Draw_AverageSignal_CollAperture.o
MAIN_RHUL_AVERSIGNAL_HALFAPERTURE = ../objects/main_Draw_AverageSignal_HalfAperture.o
MAIN_RHUL_AVERSIGNAL_ASYMMAPERTURE = ../objects/main_Draw_AverageSignal_AsymmAperture.o
MAIN_RHUL_AVERSIGNAL_PRESSURE = ../objects/main_Draw_AverageSignal_CollAperture_VacuumPressures.o
MAIN_RHUL_VOLTAGENORM = ../objects/main_Draw_AverageSignal_Voltage.o
MAIN_RHUL_AVERNOISE = ../objects/main_Draw_AverageNoise_Voltage.o
USEFUL = ../objects/UsefulFunctions.o
STYLE = ../objects/Style.o
STYLEUTILS = ../objects/Style_Utils.o

INCLUDES = -I ../include/ -I ../test/include/ -I ../../General/include/

OBJECTS = $(CONFIG) $(USEFUL) $(STYLE) $(STYLEUTILS) 

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

#../objects/test.o : ../src/ConfigReaderAnalysis.cc ../src/ConfigReaderHisto.cc ../src/ConfigReader.cc ../test/test.cc ../src/EventLooper.cc $(GTEST_HEADERS)
#	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c ../test/test.cc $(INCLUDES) -o $@

../objects/test.o : ../test/src/test.cc $(GTEST_HEADERS)
	@$(CXX) $(CPPFLAGS) -c ../test/src/test.cc $(INCLUDES) -o $@ $(CXXFLAGS)

Test : ../objects/test.o ../objects/gtest.a $(OBJECTS)
	@$(CXX) $(CPPFLAGS) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

# Building the graviton code

../BDSIM_ParticleCount_Plotter : $(MAIN_PRIMESECOND) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../BDSIM_xy_PrimariesSecondaries_Plotter : $(MAIN_XY_PRIMESECOND) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../RHULcherenkov_AverageSignal_Aperture_Plotter : $(MAIN_RHUL_AVERSIGNAL) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../RHULcherenkov_AverageSignal_HalfAperture_Plotter : $(MAIN_RHUL_AVERSIGNAL_HALFAPERTURE) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../RHULcherenkov_AverageSignal_AsymmAperture_Plotter : $(MAIN_RHUL_AVERSIGNAL_ASYMMAPERTURE) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../RHULcherenkov_AverageSignal_Aperture_VacuumPressures_Plotter : $(MAIN_RHUL_AVERSIGNAL_PRESSURE) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../RHULcherenkov_AverageSignal_VoltageNormalization_Plotter : $(MAIN_RHUL_VOLTAGENORM) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../RHULcherenkov_AverageNoise_Voltage_Plotter : $(MAIN_RHUL_AVERNOISE) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../objects/main_BDSIM_ParticleCount.o : ../src/main_BDSIM_ParticleCount.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_BDSIM_xy_PrimariesSecondaries.o : ../src/main_BDSIM_xy_PrimariesSecondaries.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_Draw_AverageSignal_CollAperture.o : ../src/main_Draw_AverageSignal_CollAperture.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_Draw_AverageSignal_HalfAperture.o : ../src/main_Draw_AverageSignal_HalfAperture.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_Draw_AverageSignal_AsymmAperture.o : ../src/main_Draw_AverageSignal_AsymmAperture.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_Draw_AverageSignal_CollAperture_VacuumPressures.o : ../src/main_Draw_AverageSignal_CollAperture_VacuumPressures.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_Draw_AverageSignal_Voltage.o : ../src/main_Draw_AverageSignal_Voltage.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/main_Draw_AverageNoise_Voltage.o : ../src/main_Draw_AverageNoise_Voltage.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/UsefulFunctions.o : ../../General/src/UsefulFunctions.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Style.o : ../../General/src/Style.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Style_Utils.o : ../../General/src/Style_Utils.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

