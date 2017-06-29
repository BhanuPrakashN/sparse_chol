CC = g++
FLAGS = -Wall -std=c++11
EXEC = smatrix.out
all:main.o smatrix.o
	$(CC) $(FLAGS) -o $(EXEC)  main.o smatrix.o 

main.o: smatrix.h main.cpp 
	$(CC) $(FLAGS) -c main.cpp 



smatrix.o: smatrix.h smatrix.cpp
	$(CC) $(FLAGS) -c smatrix.cpp



clean:
	rm *.o *.out


