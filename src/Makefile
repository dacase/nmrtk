############################################################################
#
# Makefile for GNU g++ under Linux
#
############################################################################


##############################################################################
##  Definitions of Pathes; Change as Needed for Local Site                  ##
##    In particular, BASEDIR, GALIBS, and GSLLIBS may need to be set        ##
##############################################################################
#
# For simplicity, the installed folder is assumed to be the folder of nmrtk package 
#
BASEDIR=$(HOME)/nmrtk
MYINCS=$(BASEDIR)/include
MYSRCS=$(BASEDIR)/src
MYOBJS=$(BASEDIR)/obj
MYLIBS=$(BASEDIR)/lib
MYBINS=$(BASEDIR)/bin

# For Genetic Algrithms Libs
#
# GAINCS=$(MYINCS)
GALIBS=/usr/lib/x86_64-linux-gnu

# For GNU Scientific Libs
# 
# GSLINCS=$(MYINCS)
# GSLLIBS=$(MYLIBS)
# GSLLIBNAME=libgsl.la

# For C++  libs
# Set up if the related C++ libs can not be found automatically 
CPPLIBS= 

# put all inc, src, and libs togather

CC=g++

INCS= $(CPPINCS) -I$(GAINCS) -I$(GSLINCS) -I$(MYINCS)
LIBS= $(CPPLIBS) -L$(GALIBS) -L$(GSLLIBS) -L$(MYLIBS)
CPPFLAGS=  -O3 -fpermissive -Wall $(INCS) $(LIBS) 
#CPPFLAGS= -g -Wall -fpermissive $(INCS) $(LIBS) 

LDFLAGS=  -lstdc++ -lga -lgsl -lgslcblas -lm

MAKEFILE= Makefile

SRCS=   multiFitMC.cpp \
        multiFitGA.cpp \
        UniqueSearch.cpp\
	nmr.cpp \
	relaxParm.cpp \
        calcRelaxParm.cpp \
 	relaxParmFit.cpp \
        relaxParmMultiFit.cpp \
        relaxParm2StateFit.cpp \
        S2MultiFit.cpp  

EXES=   multiFitMC \
	multiFitGA \
        UniqueSearch \
        calcRelaxParm \
	relaxParmFit \
        relaxParmMultiFit \
        relaxParm2StateFit \
	S2MultiFit     

OBJS=   nmr.o \
        relaxParm.o 
        

EXOBJS= multiFitMC.o \
	multiFitGA.o \
        calcRelaxParm.o \
	relaxParmFit.o \
	relaxParmMultiFit.o \
	relaxParm2StateFit.o \
	S2MultiFit.o 

.SUFFIXES : .o .cpp 

.cpp.o:
	$(CC) $(CPPFLAGS) -c $<

all: $(EXOBJS) $(EXES) install

bin: $(EXES)
install: 
	cp $(EXES) $(MYBINS)/
$(EXOBJS):
	$(CC) $(CPPFLAGS) -c $*.cpp

#$(EXES): $(OBJS)
#	$(CC) $(CPPFLAGS) -o $@ $@.o $(OBJS) $(LDFLAGS)
multiFitMC: multiFitMC.o
	$(CC) $(CPPFLAGS) -o multiFitMC multiFitMC.o $(LDFLAGS)
multiFitGA: multiFitGA.o
	$(CC) $(CPPFLAGS) -o multiFitGA multiFitGA.o $(LDFLAGS)
UniqueSearch: UniqueSearch.o
	$(CC) $(CPPFLAGS) -o UniqueSearch UniqueSearch.o $(LDFLAGS)    
calcRelaxParm: calcRelaxParm.o relaxParm.o nmr.o
	$(CC) $(CPPFLAGS) -o calcRelaxParm calcRelaxParm.o relaxParm.o nmr.o $(LDFLAGS)
relaxParmFit: relaxParmFit.o relaxParm.o nmr.o
	$(CC) $(CPPFLAGS) -o relaxParmFit relaxParmFit.o relaxParm.o nmr.o $(LDFLAGS)
relaxParmMultiFit: relaxParmMultiFit.o relaxParm.o nmr.o
	$(CC) $(CPPFLAGS) -o relaxParmMultiFit relaxParmMultiFit.o relaxParm.o nmr.o $(LDFLAGS)
relaxParm2StateFit: relaxParm2StateFit.o relaxParm.o nmr.o
	$(CC) $(CPPFLAGS) -o relaxParm2StateFit relaxParm2StateFit.o relaxParm.o nmr.o $(LDFLAGS)
S2MultiFit: S2MultiFit.o
	$(CC) $(CPPFLAGS) -o S2MultiFit S2MultiFit.o $(LDFLAGS)

clean:
	rm -f *.o $(EXES)

#include depend

dep: $(SRCS)
	makedepend -f depend -Y  $(SRCS)

multiFitMC.o: multiFitMC.cpp
multiFitGA.o: multiFitGA.cpp
UniqueSearch.o: UniqueSearch.cpp
calcRelaxParm.o:  calcRelaxParm.cpp relaxParm.cpp nmr.cpp
relaxParmFit.o:  relaxParmFit.cpp relaxParm.cpp nmr.cpp 
relaxParmMultiFit.o:  relaxParmMultiFit.cpp relaxParm.cpp nmr.cpp
relaxParm2StateFit.o:  relaxParm2StateFit.cpp relaxParm.cpp nmr.cpp
S2MultiFit.o:  S2MultiFit.cpp
