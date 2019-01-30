CC=gcc
CFLAGS= -Wall -Wno-float-equal -Wno-implicit-int -DAXI=1 -DMATERIALTYPE=1 -fopenmp -g 
OBJ = $(wildcard ./../../src/*.c) $(wildcard *.c) 
LIBS = -lm  -ltriangle -lmeschach  
IDIR = -I/home/stephen/Documents/software/Meschach -I/home/stephen/Documents/software/voronoi-master -I./../../headers
LIBSDIR = -L/home/stephen/Documents/software/triangle/ -L/home/stephen/Documents/software/Meschach/ 
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

simple: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(IDIR) $(LIBSDIR) $(LIBS)
