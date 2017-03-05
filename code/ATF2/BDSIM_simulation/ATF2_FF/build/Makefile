# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include -g -Wall -Wextra -pthread -O2 -std=c++1y $(shell root-config --cflags) -Wfatal-errors 

# Flags passed to the C++ compiler.
CXXFLAGS +=  $(shell root-config --glibs) -l TreePlayer

# Executables produced by this Makefile
RUN = ../CollimatorProcessAnaPlotter


# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = #Test

# Useful variables
COLLIMATORPROCESS = ../objects/Collimator_Process_Analysis.o
USEFUL = ../objects/UsefulFunctions.o
STYLE = ../objects/Style.o
STYLEUTILS = ../objects/Style_Utils.o

INCLUDES = -I ../include/ -I ../test/include/ -I ../../../../General/include/

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

../CollimatorProcessAnaPlotter : $(COLLIMATORPROCESS) $(OBJECTS)
	@$(CXX) -lpthread $^ -o $@ $(INCLUDES) $(CXXFLAGS)

../objects/ConfigReader.o : ../src/ConfigReader.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/ConfigReaderAnalysis.o : ../src/ConfigReaderAnalysis.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Collimator_Process_Analysis.o : ../src/Collimator_Process_Analysis.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/UsefulFunctions.o : ../../../../General/src/UsefulFunctions.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Style.o : ../../../../General/src/Style.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@

../objects/Style_Utils.o : ../../../../General/src/Style_Utils.cc
	@$(CXX) $(CPPFLAGS)  $(INCLUDES) -c $< -o $@
