# for gfortran
#
BINDIR="../bin"

FLAGSALWAYS = -Wall -ffree-form -DDOUBLE
FLAGSOPT= -O3 
FLAGSDEBUG= -g -fbounds-check -O0

F90=gfortran -DGFORTRAN -fopenmp
F90=mpif90 -DMPI 

#LIBS=-I/usr/local/gfortran/lib/
#-L/usr/local/openmpi-1.6.3-mx-intel12/lib/
#INC=-I/usr/local/lib/gcc/i686-apple-darwin11/4.2.1/finclude/
#-I/usr/local/include/
# choose debugging or optimization
FLAGS= ${FLAGSMACHINE} ${FLAGSALWAYS} ${INC} ${LIBS} ${FLAGSOPT}  #  change the last to ${FLAGSOPT} or ${FLAGSDEBUG}


MY_DIR=`basename ${PWD}`
#------------------------------------

%.o: %.f90
	$(F90) $(FLAGS) -c $< -o $@

%.o: %.F90
	$(F90) $(FLAGS) -c $< -o $@

all:	2pcf 3pcf 3pcf_eq 3pcf_all

2pcf: 2pcf.o kdtree2.o extension.o jackknife.o
	${F90} ${FLAGS} -o 2pcf 2pcf.o kdtree2.o extension.o jackknife.o ; \
	mv $@ $(BINDIR)

2pcf.o: kdtree2.o extension.o jackknife.o

2pcf_mix: 2pcf_mix.o kdtree2.o extension.o
	${F90} ${FLAGS} -o 2pcf_mix 2pcf_mix.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

2pcf_mix.o: kdtree2.o extension.o

3pcf: 3pcf.o kdtree2.o extension.o
	${F90} ${FLAGS} -o 3pcf 3pcf.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

3pcf.o: kdtree2.o

3pcf_eq: 3pcf_eq.o kdtree2.o extension.o
	${F90} ${FLAGS} -o 3pcf_eq 3pcf_eq.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

3pcf_eq.o: kdtree2.o

3pcf_eq2: 3pcf_eq.o kdtree2.o extension.o
	${F90} ${FLAGS} -o 3pcf_eq2 3pcf_eq2.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

3pcf_eq2.o: kdtree2.o

graph: graph.o kdtree2.o extension.o
	${F90} ${FLAGS} -o graph graph.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

graph.o: kdtree2.o

3pcf_all: 3pcf_all.o kdtree2.o extension.o
	${F90} ${FLAGS} -o 3pcf_all 3pcf_all.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

3pcf_all.o: kdtree2.o

add_fiber: add_fiber.o kdtree2.o extension.o
	${F90} ${FLAGS} -o add_fiber add_fiber.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

add_fiber.o: kdtree2.o

velocity: velocity.o kdtree2.o extension.o
	${F90} ${FLAGS} -o velocity velocity.o kdtree2.o extension.o ; \
	mv $@ $(BINDIR)

velocity.o: kdtree2.o

clean:
	/bin/rm -f *.o *.mod
