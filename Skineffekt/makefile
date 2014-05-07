CC=gcc
CFLAGS=-c -Wall -pedantic -ansi
LDFLAGS=-lm
EXECUTABLE=numerik_2

all: build

build: numerik_bespin_deutsch_2.o
	$(CC) numerik_bespin_deutsch_2.o -o $(EXECUTABLE) $(LDFLAGS)

numerik_deutsch_2.o: numerik_bespin_deutsch_2.c
	$(CC) $(CFLAGS) numerik_bespin_deutsch_2.c

clean:
	rm -rf *.o $(EXECUTABLE)
