CC=g++

FFLAGS= -O3 -g -lm 

all: 3dew.cpp
	$(CC) $(FFLAGS) 3dew.cpp -o 3dew

verify:verify.cpp
	$(CC) -O3 -g verify.cpp -o verify -lm

clean:
	rm -f 3dew
	rm -f verify
