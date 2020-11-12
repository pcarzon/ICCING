CC = g++
DEBUG = -g
CFLAGS = -std=c++17 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
MALLOC_CHECK = 2

CPP_FILES = io.cpp sorter.cpp main.cpp
HEADER_FILES = io.h event.h sorter.h
OBJECT_FILES = io.o sorter.o main.o

io.o : io.h event.h io.cpp
	$(CC) $(CFLAGS) io.cpp

sorter.o : sorter.h event.h sorter.cpp
		$(CC) $(CFLAGS) sorter.cpp

main.o : $(HEADER_FILES) main.cpp
	$(CC) $(CFLAGS) main.cpp

analysis :  $(OBJECT_FILES)
	$(CC) $(LFLAGS) $(OBJECT_FILES) -o analysis

clean :
	rm -f $(OBJECT_FILES) analysis
	echo Clean done
