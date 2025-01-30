# CWR 2024

This repository contains the course work from the 'Computergestütztes Wissenschaftliches Rechnen' lecture by Prof. S. Manmana held in the summer semester 2024.
Each numbered directory contains the solution to an exercise, along with a plotting file and makefile. The MyNumerics directory contains the self written numerical methods developed during the semester.

## The makefiles

To execute the C code use 
```
make clean
make all
```
`Note`: under windows there might be issues with the LGSL library and this might not work, depending on the setup. Everything was done on Linux.
To generate the plots use
```
make plot
```
`Note`: This uses the comand python file.py. Depending on the setup, one might have to use python3 file.py.
