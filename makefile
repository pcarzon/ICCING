CC = g++
DEBUG = -g
CFLAGS = -std=gnu++11 -Wall -c -U__STRICT_ANSI__ $(DEBUG)
LFLAGS = -Wall $(DEBUG)
MALLOC_CHECK = 2

CPP_FILES = ecc.cpp eos.cpp event.cpp io.cpp functions.cpp main.cpp splitting.cpp correlation.cpp
HEADER_FILES = ecc.h eos.h event.h io.h functions.h splitting.h global.h correlation.h
OBJECT_FILES = ecc.o eos.o event.o io.o functions.o main.o splitting.o correlation.o

event.o : event.h global.h event.cpp
	$(CC) $(CFLAGS) event.cpp

ecc.o : ecc.h event.h ecc.cpp
	$(CC) $(CFLAGS) ecc.cpp

eos.o : eos.h event.h eos.cpp
	$(CC) $(CFLAGS) eos.cpp

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
