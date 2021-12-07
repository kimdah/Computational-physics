# use this by doing this in terminal:
# make
# make run (or make compile, make clean etc.)

CPPflags = g++ -std=c++14 -O2
LIB =  -larmadillo
all: compile link clean

compile:
	${CPPflags} -c main.cpp src/Crank.cpp -Xpreprocessor -fopenmp -I include/ # ${LIB}

link:
	${CPPflags} -o main.exe main.o Crank.o -larmadillo -fopenmp #${LIB}

run:
	./main.exe

clean:
	rm -f *.o *~

plot:
	#Animation format is: animation.py 'filename' 'number of snapshots' 'snapshot times'
	python3 plotfiles/animation.py 'Problem_7.10_outputCube_slits_2.dat' 0 
	# python3 plotfiles/animation.py 'Problem_8.00_outputCube_slits_2.dat' 0.00 0.5 1.00
	# python3 ./plotfiles/animation.py 'prob9_slits_1.dat'
	# python3 ./plotfiles/animation.py 'prob9_slits_2.dat'
	# python3 ./plotfiles/animation.py 'prob9_slits_3.dat'


# CPPflags = c++ #-O3 -std=c++11 -lomp
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





#
#
#
# # --- Alternative: wildcard including all ---#
# all_alternative: compile_alternative link_alternative clean
#
# compile_alternative:
# 	${CPPflags} -c ${PROG}.cpp $(wildcard *.cpp) ${LIB}
#
# link_alternative:
# 	${CPPflags} -o ${PROG}.out $(wildcard *.o) ${LIB}
#
#
# # --- Compact: more simple ---#
# all_compact: compile_compact
#
# compile_compact:
# 	c++ -o main.out main.cpp functions.cpp
#


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
