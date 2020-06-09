CC = icpc
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o funzioneBase.o variational.o
	$(CC) random.o funzioneBase.o variational.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioneBase.o : funzioneBase.cpp funzioneBase.h
	$(CC) -c funzioneBase.cpp -o funzioneBase.o $(CFLAGS)
variational.o : variational.cpp variational.h
	$(CC) -c variational.cpp -o variational.o $(CFLAGS)
clean :
	rm *.o main.exe *.out
