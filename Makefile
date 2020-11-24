fannot:	fannot.cpp
	g++ -I"eigen.3.3.4/" -Wall -std=c++11 -g -O3 fannot.cpp -o fannot -lm	
clean:
	rm fannot


