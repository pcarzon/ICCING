CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
MALLOC_CHECK = 2

CPP_FILES = ecc.cpp eos.cpp event.cpp io.cpp main.cpp probdist.cpp splitting.cpp
HEADER_FILES = ecc.h eos.h event.h io.h probdist.h splitting.h
OBJECT_FILES = ecc.o eos.o event.o io.o main.o probdist.o splitting.o

event.o : event.h
				$(CC) $(CCFLAGS) ecc.cpp

ecc.o : ecc.h event.h
				$(CC) $(CCFLAGS) ecc.cpp

eos.o : eos.h event.h
				$(CC) $(CCFLAGS) eos.cpp

io.o : io.h event.h
				$(CC) $(CCFLAGS) io.cpp

probdist.o : probdist.h event.h
	 						$(CC) $(CCFLAGS) probdist.cpp

splitting.o : splitting.h event.h
							$(CC) $(CCFLAGS) splitting.cpp

main.o : $(HEADER_FILES) main.cpp
							$(CC) $(CCFLAGS) main.cpp

iccing :  $(CPP_FILES) $(HEADER_FILES) $(CC) $(CCFLAGS) $(OBJECT_FILES) -o iccing
