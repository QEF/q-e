
include ../make.sys

all: moduledep.x sizeof.x

moduledep.x : moduledep.c
	$(CC) $(CCFLAGS) -o moduledep.x  moduledep.c 

sizeof.x : sizeof.c
	$(CC) $(CCFLAGS) -o sizeof.x  sizeof.c

clean_:
	- rm -f moduledep.x sizeof.x clib.a *.o *.mod *.i core*

