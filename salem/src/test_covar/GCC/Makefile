CC = g++
CFLAGS = -O2 -Wall -I../../../lib/zmatlib/include

LIBS = -L../../../lib/zmatlib/ -lzmat
SRC = varcovar.cpp
OBJ = varcovar.o
MAIN = varcovar

all:    clean $(MAIN)

$(MAIN): $(OBJ) 
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJ) $(LFLAGS) $(LIBS)

$(OBJ):
	$(CC) $(CFLAGS) -c $(SRC) -o $(OBJ)

clean:
	-rm -f $(OBJ) $(MAIN)

