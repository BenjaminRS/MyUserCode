ARCH         := $(shell $(ROOTSYS)/bin/root-config --arch)
ROOTCFLAGS    = $(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTLIBS      = $(shell ${ROOTSYS}/bin/root-config --libs)
ROOTGLIBS     = $(shell ${ROOTSYS}/bin/root-config --glibs)
CERNLIB       = /cern/pro/lib
CXX = g++
LD  = g++
CXXFLAGS      = -g -O2 -Wall -fPIC -pg
LDFLAGS       = -g
SOFLAGS       = -shared
CXX           += -I${PWD}/.
CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS)
NGLIBS         = $(ROOTGLIBS)
NGLIBS        += -lMinuit -lHtml -lTreePlayer -lEG
GLIBS          = $(filter-out -lNew, $(NGLIBS))
CLIBS       = -L$(CERNLIB) -lpacklib -lmathlib -lkernlib -lpawlib -lg2c -lnsl

.SUFFIXES: .cc,.C
# ================================================================================
run:	test
test: MyObjs.o test.o ${PWD}/MyObjs.so
	$(LD) $(LDFLAGS) -o $@ $^ $(GLIBS)
${PWD}/MyObjs.so:  MyObjs.o Dict.o
	$(LD) $(SOFLAGS) -o $@ $^ $(GLIBS)
MyObjs.o:         MyObjs.C MyObjs.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
Dict.C:	MyObjs.h
	rootcint -f Dict.C -c MyObjs.h MyObjsLinkDef.h
test.o: test.cc
	$(CXX) $(CXXFLAGS) -c test.cc
clean:
	rm -rf  *.o Dict.C MyObjs.so test Dict.h 

