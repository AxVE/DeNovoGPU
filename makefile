CC = g++
CFLAGS = --std=c++14 --pedantic -pthread -Wall -I Includes/
LFLAGS = -lOpenCL

#Directories
SRCDIR = Src
OBJDIR = Obj
BINDIR = Bin
INCDIR = Includes

# Build commands

all: $(SRCDIR)/main.cpp $(OBJDIR)/workerCL.o $(OBJDIR)/log.o $(OBJDIR)/readsTools.o $(OBJDIR)/reads.o
	@mkdir -p Bin
	$(CC) $(CFLAGS) $(LFLAGS) $^ -o Bin/denovoGPU

$(OBJDIR)/workerCL.o: $(SRCDIR)/workerCL.cpp
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LFLAGS) $^ -c -o $@

$(OBJDIR)/log.o: $(SRCDIR)/log.cpp
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LFLAGS) $^ -c -o $@

$(OBJDIR)/reads.o: $(SRCDIR)/reads.cpp
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LFLAGS) $^ -c -o $@

$(OBJDIR)/readsTools.o: $(SRCDIR)/readsTools.cpp
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(LFLAGS) $^ -c -o $@

.PHONY: clean javel rebuild

clean:
	rm Obj/*

javel: clean
	rm Bin/*

rebuild: javel all
