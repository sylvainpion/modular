#
# Choose your compiler.
#
CC = gcc
# CC = g++
# CC = cc
# CC = CC

#
# Compiler flags
#
# For gcc/g++ 
CFLAGS = -O9 -Wall -Winline
#CFLAGS = -g -Wall -Winline -DDEBUG
#
# For CC
# CFLAGS = -pic -O5 +w
#
# For cc
# CFLAGS = -O5 -Dinline= -Dconst=

# For ANSI check, with gcc.
#CFLAGS = -fpic -g -W -Wall -pedantic -ansi -Winline -Dinline=__inline__

#
# Link options.
#
LD = $(CC)
#LDFLAGS = -L`pwd` -R`pwd`
LDFLAGS = -L. -lm

#
# For thread support
#CFLAGS += -DMOD_THREADS=2 -D_REENTRANT -D_MIT_POSIX_THREADS -DDEBUG
#LDFLAGS += -lpthread

###############################################
# Nothing has to be modified below this line. #
###############################################

.SUFFIXES:.c
.c.o:
	$(CC) $(INCLUDE) $(CFLAGS) -c $<

#all: libMOD.a libMOD.so time_mod.dynamic time_mod time_sign_rec
all: time_mod time_sign_rec

MOD_init.o: MOD_primes.h

gen_primes: gen_primes.o MOD_is_prime.o
	$(LD) $(LDFLAGS) -o $@ $?

MOD_primes.h: gen_primes
	./gen_primes 200 > MOD_primes.h

libMOD.so: MOD_det.o MOD_init.o MOD_is_prime.o MOD_reconstruct.o
	$(LD) $(LDFLAGS) -shared -o $@ $?
#	$(LD) $(LDFLAGS) -G -o $@ $?

libMOD.a: MOD_det.o MOD_init.o MOD_is_prime.o MOD_reconstruct.o
	ar rv $@ $+
	ranlib $@

time_mod.dynamic: time_mod.o libMOD.so
	$(LD) $(LDFLAGS) -o $@ time_mod.o -lMOD

time_mod: time_mod.o libMOD.a
	$(LD) $(LDFLAGS) -o $@ $+

time_sign_rec: time_sign_rec.o libMOD.a
	$(LD) $(LDFLAGS) -o $@ $+

clean:
	rm -f *.o gen_primes MOD_primes.h time_mod time_mod.dynamic time_mod.static time_sign_rec libMOD.a libMOD.so core
