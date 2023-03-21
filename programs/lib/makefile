SOURCES= mallocutil.c minmax.c sbess_ratio.c clib.c clear.c splint.c \
	spline.c hunt.c integral.c

OBJECTS= $(SOURCES:.c=.o)

#$(OBJECTS): $($$@:%.o=%.c)
#	$(COMPILE.c) $? -o $@
#libsmall.a: $(OBJECTS) 
#	ar cr ./libsmall.a $(OBJECTS)
#	ranlib ./libsmall.a
#	rm -f *.o
libsmall.a:  mallocutil.o minmax.o sbess_ratio.o clib.o clear.o splint.o \
        spline.o hunt.o integral.o
	ar cr ./libsmall.a mallocutil.o minmax.o sbess_ratio.o clib.o clear.o splint.o \
        spline.o hunt.o integral.o
	ranlib ./libsmall.a
	rm -f *.o
