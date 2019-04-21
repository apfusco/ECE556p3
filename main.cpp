// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv)
{

    ////  PARAMETERS  //////////////////////////////////////////////////////
    const int RUN_DURATION_IN_SEC = 240;
    ////////////////////////////////////////////////////////////////////////


  
    ////  COMMAND LINE ARGUMENTS  //////////////////////////////////////////
 	if(argc!=5){
 		printf("Usage : ./ROUTE.exe -d=0 -n=0 <input_benchmark_name> <output_file_name> \n");
 		return 1;
 	}

 	int status;

	// Get a 1 or 0 for net decomposition
	int d = atoi(&argv[1][3]);

	// Get a 1 or 0 for RRR/Maze Routing
	int n = atoi(&argv[2][3]);

	// Other CLAs
	char *inputFileName = argv[3];
 	char *outputFileName = argv[4];
	////////////////////////////////////////////////////////////////////////


	
	////  INITIAL SETUP  ///////////////////////////////////////////////////
	// create a new routing instance
 	routingInst *rst = new routingInst;  // Holds the best instance seen so far
	
 	// read benchmark
 	status = readBenchmark(inputFileName, rst);
 	if(status==0){
 		printf("ERROR: reading input file \n");
 		return 1;
 	}
	////////////////////////////////////////////////////////////////////////

	

	////  NET DECOMPOSITION (if necessary)  ////////////////////////////////
	if(d == 1){
	  status = netDecomp(rst);
	  if(status==0){
	    printf("ERROR: decomposing nets\n");
	    release(rst);
	    return 1;
	  }
	}
	////////////////////////////////////////////////////////////////////////

	
	
 	////  CREATE INITIAL ROUTING SOLUTION  /////////////////////////////////
 	status = solveRoutingBasic(rst);
 	if(status==0){
 		printf("ERROR: running routing \n");
 		release(rst);
 		return 1;
 	}	
	////////////////////////////////////////////////////////////////////////



	////  RIPUP AND RE-ROUTE (if necessary)  ///////////////////////////////
	if (n == 1) {
	  // Holds the time
	  clock_t start, now;
	  double time_taken;
	  
	  // Get overflow of the initial solution
	  int lowestOverFlow = getOverflow(rst);
	  int currOverFlow;

	  // routing instance to route next
	  routingInst *curr_rst;

	  // Get current time
	  start = clock();
	  
	  // Main RRR Loop (terminates after run duration is reached)
	  do {
		  // Initialize this routing instance
		  curr_rst = new routingInst;
		  
		  // read benchmark for curr_rst
	      status = readBenchmark(inputFileName, curr_rst);
   	      if(status==0){
 	  	    printf("ERROR: reading input file for curr_rst \n");
			release(curr_rst);
			release(rst);
 	        return 1;
 	      }
		  
		  // Make order array
		  int *order = makeOrderArray(rst);
		  
		  // Solve routing for curr_rst
		  status = solveRoutingRand(order, curr_rst);
		  free(order);
 	      if(status==0){
 		    printf("ERROR: running routing for curr_rst \n");
 		    release(rst);
			release(curr_rst);
 		    return 1;
 	      }
		  
		  // Check if curr_rst is the best solution we've seen so far
		  currOverFlow = getOverflow(curr_rst);
		  
		  if(currOverFlow < lowestOverFlow) {
			  // This is the best solution seen so far
			  lowestOverFlow = currOverFlow;
			  // Replace rst with curr_rst
			  release(rst);
			  rst = curr_rst;
		  } else {
			  // This solution is bad so destroy it
			  release(curr_rst);
		  }
			
	      // Find elapsed time since start
	      now = clock() - start;
	      time_taken = ((double)now)/CLOCKS_PER_SEC; // in seconds
	    
	  } while (time_taken < RUN_DURATION_IN_SEC);
	}
	////////////////////////////////////////////////////////////////////////
	
	
	
 	////  WRITE THE RESULT  ////////////////////////////////////////////////
 	status = writeOutput(outputFileName, rst);
 	if(status==0){
 		printf("ERROR: writing the result \n");
 		release(rst);
 		return 1;
 	}
	////////////////////////////////////////////////////////////////////////

	

	////  END  /////////////////////////////////////////////////////////////
 	release(rst);
 	printf("\nDONE!\n");	
 	return 0;
	////////////////////////////////////////////////////////////////////////
}
