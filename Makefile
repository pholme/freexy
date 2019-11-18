SRC = .
CFLAGS = -W -Wall -Ofast -march=native
LDFLAGS = -lm
CC = gcc

OBJ1 = o/freexy.o o/pcg_rnd.o o/updates.o o/measures.o
OBJ2 = o/freerxy.o o/pcg_rnd.o o/updates_freer.o o/measures_freer.o 

all : freexy freerxy

freexy: $(OBJ1)
	$(CC) -o $@ $^ $(LDFLAGS)

freerxy: $(OBJ2)
	$(CC) -o $@ $^ $(LDFLAGS)

o/freexy.o : $(SRC)/freexy.c $(SRC)/freexy.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/freexy.c -o $@

o/freerxy.o : $(SRC)/freexy.c $(SRC)/freexy.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -DFREER -c $(SRC)/freexy.c -o $@

o/updates.o : $(SRC)/updates.c $(SRC)/freexy.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/updates.c -o $@

o/updates_freer.o : $(SRC)/updates.c $(SRC)/freexy.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -DFREER -c $(SRC)/updates.c -o $@

o/measures.o : $(SRC)/measures.c $(SRC)/freexy.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/measures.c -o $@

o/measures_freer.o : $(SRC)/measures.c $(SRC)/freexy.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -DFREER -c $(SRC)/measures.c -o $@

o/pcg_rnd.o : $(SRC)/pcg_rnd.c $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/pcg_rnd.c -o $@

