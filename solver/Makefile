#!/usr/bin/make

CFLAGS=-O3 -Wall -std=c90
CC=gcc

OBJS = \
	daux.o ddaskr.o dlinpk.o

# target dir for install
DESTDIR=/usr/local

#  Static library target
#

all: libcdaskr.a

libcdaskr.a:  $(OBJS)
	ar r $@ $(OBJS)

%.o: %.c
	$(CC) -c -o $@ $<  $(CFLAGS)

install: libcdaskr.a
	cp libcdaskrk.a ${DESTDIR}/lib
	chmod 644 ${DESTDIR}/lib/libcdaskr.a

clean:
	  rm -f *.o libcdaskr.a *~ #*#

remove: libcdaskr.a
	rm -f ${DESTDIR}/lib/libcdaskr.a
