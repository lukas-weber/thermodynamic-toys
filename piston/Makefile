CC=gcc
CFLAGS=-O3 -ffast-math

libpiston.so : piston.c
	$(CC) $(CFLAGS) -fPIC -shared -Wl,-soname,libpiston.so -o libpiston.so piston.c
