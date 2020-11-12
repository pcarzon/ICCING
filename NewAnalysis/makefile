CC = g++
DEBUG = -g
CFLAGS = -std=c++17 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
MALLOC_CHECK = 2

CPP_FILES = io.cpp sorter.cpp calculator.cpp main.cpp
HEADER_FILES = io.h event.h sorter.h calculator.h
OBJECT_FILES = io.o sorter.o calculator.o main.o

io.o : io.h event.h io.cpp
	$(CC) $(CFLAGS) io.cpp

sorter.o : sorter.h event.h sorter.cpp
	$(CC) $(CFLAGS) sorter.cpp

calculator.o : calculator.h event.h calculator.cpp
	$(CC) $(CFLAGS) calculator.cpp

main.o : $(HEADER_FILES) main.cpp
	$(CC) $(CFLAGS) main.cpp

analysis :  $(OBJECT_FILES)
	$(CC) $(LFLAGS) $(OBJECT_FILES) -o analysis

clean :
	rm -f $(OBJECT_FILES) analysis
	echo Clean done
