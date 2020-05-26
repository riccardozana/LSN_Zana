CC = icpc
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o GA.o
	$(CC) random.o GA.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
GA.o : GA.cpp GA.h
	$(CC) -c GA.cpp -o GA.o $(CFLAGS)
clean :
	rm *.o main.exe *.out
