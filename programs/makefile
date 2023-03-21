HOME=/home/watada/x86_64/src/tsunami/programs
CPPFLAGS = -I$(HOME)/lib
LDFLAGS =  -L$(HOME)/lib
test_read_solid: test_read_solid.c read_solid.o model.o ray_sub.o \
calc_a.earth.o read_atmos.o search.o 
	$(LINK.c) test_read_solid.c read_solid.o model.o ray_sub.o \
	calc_a.earth.o read_atmos.o search.o \
	-lm -lsmall -o test_read_solid
solid_rk: solid_rk.c rk.o read_solid.o read_atmos.o calc_a.o model.o \
solid_sub.o eigen.o 
	$(LINK.c) -DDBG solid_rk.c rk.o read_solid.o read_atmos.o \
	calc_a.o model.o solid_sub.o eigen.o \
	-lm -lsmall -o solid_rk
clean:
	rm *.o solid_rk test_read_solid
