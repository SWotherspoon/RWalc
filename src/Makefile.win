all: RWalc.so

RWalc.so: RWalc.cpp
	Rscript --vanilla -e "TMB::compile('RWalc.cpp')"

clean:
	rm -rf *.o
