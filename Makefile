all: ec

#OPT=-O0 -g
OPT=-Ofast

CFLAGS=$(OPT) -Wall -std=gnu99 -D_XOPEN_SOURCE=500
LINK=-pthread -lm

ec.o: ec.c rng.h ecs.h sys.h ec_run.h
	$(CC) $(CFLAGS) -c $<

# shut up make
flt0.yy:;

flt0.yy.c flt0.yy.h: flt0.yy flt0.h
	flex -7 -Pflt0 -oflt0.yy.c --header-file=flt0.yy.h $<

flt0.yy.o: flt0.yy.c flt0.h
	$(CC) $(CFLAGS) -c $<

ec_run.o: ec_run.c rng.h sys.h ec_run.h cpu.h ecs.h flt0.yy.h flt0.h
	$(CC) $(CFLAGS) -c $<

ecs.o: ecs.c ecs.h sys.h
	$(CC) $(CFLAGS) -c $<

ec: ec.o ec_run.o ecs.o flt0.yy.o
	$(CC) $(LINK) $^ -o $@

clean:
	rm -rf *.o *.yy.c *.yy.h ec
