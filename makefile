CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
MALLOC_CHECK = 2

CPP_FILES = ecc.cpp eos.cpp event.cpp io.cpp main.cpp probdist.cpp splitting.cpp
HEADER_FILES = ecc.h eos.h event.h io.h probdist.h splitting.h
OBJECT_FILES = ecc.o eos.o event.o io.o main.o probdist.o splitting.o

event.o : event.h event.cpp
	$(CC) $(CFLAGS) event.cpp

ecc.o : ecc.h event.h ecc.cpp
	$(CC) $(CFLAGS) ecc.cpp

eos.o : eos.h event.h eos.cpp
	$(CC) $(CFLAGS) eos.cpp

io.o : io.h event.h io.cpp
	$(CC) $(CFLAGS) io.cpp

probdist.o : probdist.h event.h probdist.cpp
	$(CC) $(CFLAGS) probdist.cpp

splitting.o : splitting.h event.h splitting.cpp
	$(CC) $(CFLAGS) splitting.cpp

main.o : $(HEADER_FILES) main.cpp
	$(CC) $(CFLAGS) main.cpp

iccing :  $(OBJECT_FILES)
	$(CC) $(LFLAGS) $(OBJECT_FILES) -o iccing

clean :
	rm -f $(OBJECT_FILES) iccing
	echo Clean done
