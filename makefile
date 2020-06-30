CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -Wall -c $(DEBUG)
CXXFLAGS = $(shell pkg-config --cflags opencv)
LIBS = $(shell pkg-config --libs opencv)
LFLAGS = -Wall $(DEBUG)
MALLOC_CHECK = 2

CPP_FILES = ecc.cpp eos.cpp event.cpp io.cpp main.cpp probdist.cpp splitting.cpp
HEADER_FILES = ecc.h eos.h event.h io.h probdist.h splitting.h
OBJECT_FILES = ecc.o eos.o event.o io.o main.o probdist.o splitting.o

event.o : event.h event.cpp
	$(CC) $(CCFLAGS) event.cpp

ecc.o : ecc.h event.h ecc.cpp
	$(CC) $(CCFLAGS) ecc.cpp

eos.o : eos.h event.h eos.cpp
	$(CC) $(CCFLAGS) eos.cpp

io.o : io.h event.h io.cpp
	$(CC) $(CCFLAGS) io.cpp

probdist.o : probdist.h event.h probdist.cpp
	$(CC) $(CCFLAGS) probdist.cpp

splitting.o : splitting.h event.h splitting.cpp
	$(CC) $(CCFLAGS) splitting.cpp

main.o : $(HEADER_FILES) main.cpp
	$(CC) $(CCFLAGS) main.cpp

iccing :  main.o
	$(CC) $(CCFLAGS) main.o -o iccing

clean :
	rm -f $(OBJECT_FILES)
	echo Clean done
