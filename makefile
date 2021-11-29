
CPPflags = g++ -std=c++11
LIB = -I /usr/local/opt/armadillo/include -DARMA_USE_LAPACK -DARMA_USE_BLAS -DARMA_DONT_USE_WRAPPER -llapack -lblas

all: compile link clean

compile:
	${CPPflags} -c main.cpp PenningTrap.cpp Particle.cpp

link:
	${CPPflags} -o main.out main.o ${LIB}

run:
	./main.out

clean:
	rm -f *.o *~

# CPPflags = c++ #-O3 -std=c++11
# LIB =  #-larmadillo # -llapack -lblas
# PROG = main
#
# all: compile link clean #run
#
# compile:
# 	${CPPflags} -c ${PROG}.cpp functions.cpp ${LIB}
#
# link:
# 	${CPPflags} -o ${PROG}.out ${PROG}.o functions.o ${LIB}
#
# run:
# 	./${PROG}.out
#
# clean:
# 	rm -f *.o *~








# --- Alternative: wildcard including all ---#
#all_alternative: compile_alternative link_alternative clean

#compile_alternative:
#	${CPPflags} -c ${PROG}.cpp $(wildcard *.cpp) ${LIB}

#link_alternative:
#	${CPPflags} -o ${PROG}.out $(wildcard *.o) ${LIB}


# --- Compact: more simple ---
#all_compact: compile_compact

#compile_compact:
#	c++ -o main.out main.cpp functions.cpp



#--- Old method ---#
# CPPflags = c++ #-O3 -std=c++11
# LIB =  #-larmadillo -llapack -lblas
# PROG = main
#
#
# ${PROG}: 						${PROG}.o functions.o
# 										${CPPflags} ${PROG}.o functions.o ${LIB} -o ${PROG}.exe
#
#
# ${PROG}.o: 					${PROG}.cpp
# 										${CPPflags} -c ${PROG}.cpp
#
# functions.o: 				functions.cpp functions.hpp
# 										${CPPflags} -c functions.cpp
#
#
# clean:
# 	rm -f *.o *~
