CC = gcc

LIB = /cs/labs/dina/dina/libs/gamb++/
MYLIB = /cs/labs/dina/dina/libs/gamb++/../DockingLib/

CCFLAGS = -I$(LIB) -I$(MYLIB) -Wall -O2
LDFLAGS = -static -L$(LIB) -L$(MYLIB) -ldockingLib -lgamb++ -lstdc++ -O2 -lm

# User defined classes and modules. (no file suffixes)

CLASSES = MolecularInterface interface

# Prepare object and source file list using pattern substitution func.
ALL  = $(CLASSES)
OBJS = $(patsubst %, %.o,  $(ALL))
SRCS = $(patsubst %, %.cc, $(ALL))

TARGET = interface

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)  -o $(TARGET)

%.o: %.cc
	$(CC) $(CCFLAGS) -c $*.cc

clean:
	/bin/rm -f *.o *~ \#* core

backup:
	zip rmsd.zip *.cc *.h Makefile* readme *.txt run

depend:
	makedepend -- $(CCFLAGS) -- $(SRCS)
# DO NOT DELETE THIS LINE -- make depend depends on it.
