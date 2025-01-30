
default:
	gcc -o output.exe exercise5.c my_numerics.c -lm

run:
	./output.exe

clean:
	rm *.o
	rm *.exe

exercise5: exercise5.o my_numerics.o
	gcc exercise5.o my_numerics.o -lm -o exercise5

exercise5.o: my_numerics.h exercise5.c
	gcc -c exercise5.c

my_numerics.o: my_numerics.h my_numerics.c
	gcc -c my_numerics.c