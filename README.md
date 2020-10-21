# RLibm

RLibm is a math library that provides correctly rounded result for all inputs. Currently, RLibm supports a number of elementary functions for bfloat16, posit16, and float representations. rlibm is generated using the technique described from this [paper](https://arxiv.org/pdf/2007.05344.pdf).

## Installation
To compile the math library, please follow the instructions below.
This compilation instruction conrates separate math library for each of the available representations.

### Prerequisite
If you want to compile the math library for posit16, you have to install SoftPosit. Please follow the instructions from the [SoftPosit GitLab](https://gitlab.com/cerlane/SoftPosit).

### Installation step
1. Clone the rlibm repository
```
git clone https://github.com/rutgers-apl/rlibm.git
```

2. Create an environment variable SOFTPOSITPATH that points to the directory of SoftPosit:
```
export SOFTPOSITPATH=<path-to-softposit-directory>
```
  
3. Build the math library
  1. If you want to build all the math libraries, simply use make rule
  ```
  cd rlibm
  make
  ```

  2. If you want to build math libraries for each representation separately, you can use these make rule
  ```
  cd rlibm
  make bfloat16mlib
  make posit16mlib
  make floatmlib
  ```
4. (Optional) Test that the math library does produce the correct value. This step requires MPFR installed.
  1. To run the correctness bench suite for all math libraries, run
  ```
  ./runLibTest.sh
  ```
  2. The correctness bench suite for math library is located in the libtest folder.
  ```
  cd libtest/bfloat16
  make
  ./runAll.sh
  ```

## USAGE
The math library will be located in the lib directory:
  * bfloat16MathLib.a : math library for bfloat16
  * floatMathLib.a : math library for float
  * posit16MathLib.a : math library for posit16.

The header files for each library is located in the include directory:
  * bfloat16_math.hpp : header for bfloat16 math library
  * float_math.h : header for float math library
  * posit16_math.h : header for posit16 math library
  
If you want to use the bfloat16 math library, you need to also include `bfloat16.hpp` which constains our custom bfloat16 class.

You can use our library in the code similar to how standard math library is used, except our function names start with "rlibm_":
```
test.cpp: 
#include "float_math.h"
int main() {
  float result = rlibm_cospi(1.5f);
  return 0;
}
```

To build the program, include the math library in the compilation command:
```
g++ test.cpp ../../lib/floatMathLib.a -lm -o test
```
Currently, RLibm uses some functions from the default math library for range reduction, such as to decompose a floating point value into the integral part and fractional part.
