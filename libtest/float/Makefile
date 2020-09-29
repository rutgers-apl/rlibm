CC = g++

funcNames = myLog2

src = $(patsubst %, %.cpp, $(funcNames))
obj = $(patsubst %, %, $(funcNames))

all: $(obj)

%: %.cpp
	$(CC) -std=c++11  -I../../helper/ -I../../include/  $^ ../../lib/floatMathLib.a -lm -lmpfr -lgmp -o $@

clean:
	rm -rf $(obj) *.dSYM
