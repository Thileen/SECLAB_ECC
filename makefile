OBJ = testbench.o seclabec.o

prog: $(OBJ)
	gcc -o prog $(OBJ)

testbench.o: testbench.c ec_domain_parameter.h seclabec.h
	gcc -c -Wall testbench.c
	
seclabec.o: seclabec.c ec_domain_parameter.h seclabec.h
	gcc -c -Wall seclabec.c 



	