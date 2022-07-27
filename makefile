all: ball_samp
ball_samp: ball_samp.o
	g++ -fopenmp ball_samp.o -o ball_samp
ball_samp.o: ball_samp.cpp
	g++ -c ball_samp.cpp -fopenmp