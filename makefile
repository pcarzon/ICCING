CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
MALLOC_CHECK = 2

CPP_FILES = ecc.cpp eos.cpp event.cpp io.cpp main.cpp probdist.cpp splitting.cpp
HEADER_FILES = ecc.h eos.h event.h io.h probdist.h splitting.h
OBJECT_FILES = ecc.o eos.o event.o io.o main.o probdist.o splitting.o

iccing :  $(CPP_FILES) $(HEADER_FILES) $(CC) $(CCFLAGS) $(OBJECT_FILES) -o iccing
