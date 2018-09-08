CC=gcc
LIB=-lz -lm

INCLUDE=

CFLAGS=$(INCLUDE) -fPIC -Wall

objects=bedgraph.o chr.o fastq.o memutil.o nuc.o seq.o seqcoord.o txtfile.o util.o wig.o

default: all

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@ $(INCLUDE)

libgenome.so: $(objects)
	$(CC) -$(LIB) -shared $(objects) -o libgenome.so


all: $(objects) libgenome.so

clean:
	rm -f $(objects) libgenome.so
