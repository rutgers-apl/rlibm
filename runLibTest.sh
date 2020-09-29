#!/bin/bash


echo -e "\033[1mBuilding math libraries \033[0m"
make -s clean
make -s


echo -e "\033[1mTesting bfloat16 math library correctness \033[0m"
cd libtest/bfloat16
make -s clean
make -s
./runAll.sh
make -s clean


echo -e "\033[1mTesting posit16 math library correctness \033[0m"
cd ../posit16
make -s clean
make -s
./runAll.sh
make -s clean


echo -e "\033[1mTesting float math library correctness \033[0m"
cd ../float
make -s clean
make -s
./runAll.sh
make -s clean


cd ../..
