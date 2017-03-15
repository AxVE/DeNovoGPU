CC = g++
CFLAGS = --std=c++14 --pedantic -pthread -Wall -I Includes/
LFLAGS = -lOpenCL

#Directories
SRCDIR = Src
OBJDIR = Obj
BINDIR = Bin
INCDIR = Includes

#Files name(don't input main)
files = log reads



all: $(SRCDIR)/main.cpp $(OBJDIR)/workerCL.o $(OBJDIR)/log.o $(OBJDIR)/readsTools.o $(OBJDIR)/reads.o
	$(CC) $(CFLAGS) $(LFLAGS) $^ -o Bin/denovoGPU

$(OBJDIR)/workerCL.o: $(SRCDIR)/workerCL.cpp
	$(CC) $(CFLAGS) $(LFLAGS) $^ -c -o $@

$(OBJDIR)/log.o: $(SRCDIR)/log.cpp
	$(CC) $(CFLAGS) $(LFLAGS) $^ -c -o $@

$(OBJDIR)/reads.o: $(SRCDIR)/reads.cpp
	$(CC) $(CFLAGS) $(LFLAGS) $^ -c -o $@

$(OBJDIR)/readsTools.o: $(SRCDIR)/readsTools.cpp
	$(CC) $(CFLAGS) $(LFLAGS) $^ -c -o $@

.PHONY: clean javel rebuild

clean:
	rm Obj/*

javel: clean
	rm Bin/*

rebuild: javel all
