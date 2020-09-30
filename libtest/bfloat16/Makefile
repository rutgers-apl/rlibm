CC = g++

funcNames = Log2 Log Log10 Sinpi Cospi Sqrt Cbrt Exp Exp2 Exp10

src = $(patsubst %, %.cpp, $(funcNames))
obj = $(patsubst %, %, $(funcNames))

all: $(obj)

%: %.cpp
	$(CC) -std=c++14 $^ ../../lib/bfloat16MathLib.a -lm -lmpfr -lgmp -I../../helper/ -I../../include/ -o $@

clean:
	rm -rf $(obj) *.dSYM
