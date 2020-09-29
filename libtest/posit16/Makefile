objects = log2 log log10 sinpi cospi sqrt exp exp2 exp10

all: $(objects)

$(objects): %: %.cpp
	g++ -std=c++14 $^ -lm -lmpfr -lgmp ../../lib/posit16MathLib.a $(SOFTPOSITPATH)/build/Linux-x86_64-GCC/softposit.a  -I$(SOFTPOSITPATH)/source/include -I../../include/ -o $@


clean:
	rm -rf $(objects) *.dSYM
