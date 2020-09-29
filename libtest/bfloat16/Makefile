CC = g++

funcNames = myLog2 myLog myLog10 mySinpi myCospi mySqrt myCbrt myExpv2 myExp2v2 myExp10v2

src = $(patsubst %, %.cpp, $(funcNames))
obj = $(patsubst %, %, $(funcNames))

all: $(obj)

%: %.cpp
	$(CC) -std=c++14 $^ ../../lib/bfloat16MathLib.a -lm -lmpfr -lgmp -I../../helper/ -I../../include/ -o $@

clean:
	rm -rf $(obj) *.dSYM
