CC = icpc
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o functions.o funzioneBase.o integralMC.o
	$(CC) functions.o random.o funzioneBase.o integralMC.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
funzioneBase.o : funzioneBase.cpp funzioneBase.h
	$(CC) -c funzioneBase.cpp -o funzioneBase.o $(CFLAGS)
integralMC.o : integralMC.cpp integralMC.h
	$(CC) -c integralMC.cpp -o integralMC.o $(CFLAGS)
clean :
	rm *.o main.exe *.out
