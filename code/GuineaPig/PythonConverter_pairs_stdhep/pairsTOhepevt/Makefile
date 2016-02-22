#
# ROOT Makefile
#

# Commands and options for compiling
CC = g++
CXXFLAGS = -O3 -Wall -Wextra -Wshadow -fPIC  -I$(ROOTSYS)/include #-I$(HOME)/root_macros
CXXFLAGSNOROOT = -O3 -Wall -Wextra -Wshadow -fPIC #-I$(HOME)/root_macros
# Options for linking
LDFLAGS = -O `root-config --libs` -lMinuit -lMathCore
LDFLAGSNOROOT = -O 
#-L$(ROOTSYS)/lib -lNew -lCore -lCint -lHist -lGraf \
#          -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix \
#          -lPhysics -lm -ldl -rdynamic

#
# Add your files and executables here

all:pairToHepEvt


pairToHepEvt: pairToHepEvt.o
	$(CC) $(LDFLAGSNOROOT) $^ -o $@

debug: pairToHepEvt.C
	$(CC) $(CXXFLAGSNOROOT) -DDEBUG -c -o debug.o $^
	$(CC) $(LDFLAGSNOROOT) debug.o -o pairToHepEvt

#
# Cleanup
clean:
	rm -f DEADJOE *~ *.o coherentToHEPEvt pairToHepEvt pairMuons pairElectrons debug
