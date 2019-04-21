ECE 556 Spring 2019
Project - Part 2 - Global Router

Authors:
Cameron Craig - Net Decomp, RRR, Net Re-ordering, Termination Condition
Alex Fusco - Maze Routing (A*), File I/O, Memory Management

Files included:
ece556.h
ece556.cpp
main.cpp
Makefile
README.md


Purpose:
This global router takes in a file which lists nets, with points on the net which need to be routed together.
This version of the router uses a simple approach initially, then runs iterations of A* to route the pins.
At the termination point the best solution is output.


Usage:
1) Enter "make" on the command line to compile ROUTE.exe
2) Enter: ./ROUTE.exe -d=? -n=? [benchmark name] [solution file name]
	-d: set to 1 to enable net decomposition
	-n: set to 1 to enable RRR, Edge Weight Calc, and Maze Routing
3) (Optional) Enter: ./556_eval [benchmark name] [solution file name] [eval mode]
   where eval_mode=0 runs a thourough evaluation of the generated route.

This will create a file with the given solution file name that contains the computed route.


Measured Performance:

	1) -d=0 -n=0
	TWL: 5465212
	TOF: 1669008
	
	2) -d=0 -n=1
	TWL: 5487152
	TOF: 1253171
	
	3) -d=1 -n=0
	TWL: 3768801
	TOF: 698129
	
	4) -d=1 -n=1
	TWL: 3768801
	TOF: 698129
