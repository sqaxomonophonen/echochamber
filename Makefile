all: ec

#OPT=-O0 -g
OPT=-Ofast

CFLAGS=$(OPT) -Wall -std=gnu99 -D_XOPEN_SOURCE=500
LINK=-pthread -lm

ec.o: ec.c rng.h ecs.h sys.h ec_run.h
	$(CC) $(CFLAGS) -c $<

# shut up make
flt0.ll:;
flt0.ll.c flt0.ll.h: flt0.ll flt0.h
	flex -7 -Pflt0 -oflt0.ll.c --header-file=flt0.ll.h $<

flt0.ll.o: flt0.ll.c flt0.h
	$(CC) $(CFLAGS) -c $<

ec_run.o: ec_run.c rng.h sys.h ec_run.h cpu.h ecs.h flt0.ll.h flt0.h
	$(CC) $(CFLAGS) -c $<

ecs.o: ecs.c ecs.h sys.h
	$(CC) $(CFLAGS) -c $<

ec: ec.o ec_run.o ecs.o flt0.ll.o
	$(CC) $(LINK) $^ -o $@

clean:
	rm -rf *.o *.ll.c *.ll.h ec
