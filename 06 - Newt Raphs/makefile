Main = ./lagrange.c
Sources = ../MyNumerics/my_numerics.c
Output = lagrange

Libs = -lm

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
	python3 plot.py