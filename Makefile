# DSP ENGI4151 COURSEWORK 19/20 (6_2_2020)
# GROUP: 12
# NAMES: EVAN SUTCLIFFE, WILL PANTON, JAMIE MCDONALD

# define the name of your source file(s)
SRCS = group12dsp.c fft.c

# define the name of the object files(s) - we can do this automatically
OBJS = $(SRCS:.c=.o)

# tell MAKE which compiler to use
CCOMP = gcc

# flags for the compiler
CFLAGS = -Wall -O2 -fstrict-aliasing

# flags for the linker - note -lm for the math library
LDFLAGS = -O2 -lm -L/usr/lib -pthread

# the name of your executable file (the target) - here we put it in the top directory
TARGET = group12dsp

# actions
all: $(OBJS)
	$(CCOMP) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CCOMP) -c -o $@ $< $(CFLAGS)

# delete all objects and target
clean:
	rm -f $(OBJS) $(TARGET)
