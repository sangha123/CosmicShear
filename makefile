.SUFFIXES: .cpp .o
CPP = g++ -g -fopenmp 
flags = -O3 -Wall 


lflags = -lm -lgsl -lgslcblas 


all: powcamb_v0

clean:
	-rm *.o
	-rm *~

.cpp.o: $?
	$(CPP) $(flags) -c $*.cpp

powcamb_v0: powcamb_v0.o powfunc.o 
	$(CPP) powcamb_v0.o powfunc.o $(lflags) -o powcamb_v0

