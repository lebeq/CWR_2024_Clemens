Main = ./exercise7.c
Sources = ../MyNumerics/my_numerics.c
Output = exercise7

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
	rm -f *.csv
	rm -f *.png
	rm -f *.exe

plot:
	python3 plots.py