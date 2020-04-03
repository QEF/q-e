# Usage
# make        # compile all binary
# make clean  # remove ALL binaries and objects

CC=ifort   # compiler to use
CFLAGS=

OBJS=create_hz_input create_pw_input2

all: $(OBJS)

$(OBJS): %: %.f90
	$(CC) $(CFLAGS) -o $@ $<
 
clean : 
	rm -f $(OBJS)
