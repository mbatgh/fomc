IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -O2 -Wunused

ODIR=obj

LIBS=-lm

_DEPS = functions.h strucs.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

SOURCES = *.c

#_OBJ=$(SOURCES:.c=.o)

_OBJ = main.o micutc.o various.o get_parameters.o read_atp.o initialize_sigeps.o initialize_morphs.o

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

fomc:	$(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ $(IDIR)/*~ fomc *.backup

install: fomc
	cp fomc ~/bin/

