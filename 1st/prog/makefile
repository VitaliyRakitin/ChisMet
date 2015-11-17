CC=gcc
CFLAGS=-c -Wall -O3 -Werror
LDFLAGS=
SOURCES=main.c k.c special.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=p

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
	rm -rf *.o

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm -rf *.o 