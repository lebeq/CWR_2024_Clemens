Main = ./pendulums.c
Sources = ../MyNumerics/my_numerics.c
Output = pendulums

Libs = -lm -lgsl -lgslcblas

Flags = -Wall -Wextra

#All-Befehl
all: compile run

compile:
	gcc $(Main) $(Sources) -o $(Output) $(Libs) $(Flags)

run:
	./$(Output)

clean:
	rm -f $(Output)
	rm -f *.txt
	rm -f *.csv

plot:
	python3 plots.py





#PENDULUMS: pendulums.o my_numerics.o
#	gcc -Wall pendulums.o my_numerics.o -lgsl -lgslcblas -lm -o PENDULUMS
#
#pendulums.o: pendulums.c ../MyNumerics/my_numerics.h
#	gcc -I.. -Wall -c pendulums.c
#
#my_numerics.o: ../MyNumerics/my_numerics.h ../MyNumerics/my_numerics.c
#	gcc -Wall -c ../MyNumerics/my_numerics.c
