CC = g++
DEBUG = -g
CFLAGS = -std=c++17 -Wall -c -fopenmp `gsl-config --cflags` `gsl-config --libs` $(DEBUG)
LFLAGS = -Wall -fopenmp `gsl-config --cflags` `gsl-config --libs` $(DEBUG)
MALLOC_CHECK = 2

CPP_FILES = ecc.cpp event.cpp io.cpp functions.cpp main.cpp splitting.cpp correlation.cpp
HEADER_FILES = ecc.h event.h io.h functions.h splitting.h global.h correlation.h
OBJECT_FILES = ecc.o event.o io.o functions.o main.o splitting.o correlation.o

ecc.o : ecc.h ecc.cpp
	$(CC) $(CFLAGS) ecc.cpp

event.o : event.h global.h ecc.h event.cpp
	$(CC) $(CFLAGS) event.cpp

correlation.o : correlation.h correlation.cpp
	$(CC) $(CFLAGS) correlation.cpp

io.o : io.h event.h correlation.h io.cpp
	$(CC) $(CFLAGS) io.cpp

functions.o : functions.h functions.cpp
		$(CC) $(CFLAGS) functions.cpp

splitting.o : splitting.h global.h functions.h correlation.h splitting.cpp
	$(CC) $(CFLAGS) splitting.cpp

main.o : $(HEADER_FILES) main.cpp
	$(CC) $(CFLAGS) main.cpp

iccing :  $(OBJECT_FILES)
	$(CC) $(LFLAGS) $(OBJECT_FILES) -o iccing

clean :
	rm -f $(OBJECT_FILES) iccing
	echo Clean done
