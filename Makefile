CC=g++
CFLAGS=-g -O2
LFLAGS=

all: hotmix gen-rand-hot gen-data likelihood parent

common.o: common.cpp common.hpp common-impl.hpp
	$(CC) $(CFLAGS) -c common.cpp

hotmix.o: hotmix.cpp common.hpp common-impl.hpp edmonds_optimum_branching.hpp rooted-tree.hpp
	$(CC) $(CFLAGS) -c hotmix.cpp

gen-data.o: gen-data.cpp common.cpp common.hpp common-impl.hpp rooted-tree.hpp
	$(CC) $(CFLAGS) -c gen-data.cpp

gen-rand-hot.o: gen-rand-hot.cpp common.cpp common.hpp common-impl.hpp rooted-tree.hpp
	$(CC) $(CFLAGS) -c gen-rand-hot.cpp

likelihood.o: likelihood.cpp common.cpp common.hpp common-impl.hpp rooted-tree.hpp
	$(CC) $(CFLAGS) -c likelihood.cpp

parent.o: parent.cpp common.cpp common.hpp common-impl.hpp rooted-tree.hpp
	$(CC) $(CFLAGS) -c parent.cpp

hotmix: hotmix.o common.o
	$(CC) $(LFLAGS) hotmix.o common.o -lmpfrcpp -lmpfr -lgmpxx -lgmp -lNHparser -lboost_program_options -o hotmix

gen-rand-hot: gen-rand-hot.o common.o
	$(CC) $(LFLAGS) gen-rand-hot.o common.o -lmpfrcpp -lmpfr -lgmpxx -lgmp -lNHparser -lboost_program_options -o gen-rand-hot

gen-data: gen-data.o common.o
	$(CC) $(LFLAGS) gen-data.o common.o -lmpfrcpp -lmpfr -lgmpxx -lgmp -lNHparser -lboost_program_options -o gen-data

likelihood: likelihood.o common.o
	$(CC) $(LFLAGS) likelihood.o common.o -lmpfrcpp -lmpfr -lgmpxx -lgmp -lNHparser -lboost_program_options -o likelihood

parent: parent.o common.o
	$(CC) $(LFLAGS) parent.o common.o -lmpfrcpp -lmpfr -lgmpxx -lgmp -lNHparser -lboost_program_options -o parent

clean:
	rm -f *.o
	rm -f hotmix gen-rand-hot gen-data likelihood parent
