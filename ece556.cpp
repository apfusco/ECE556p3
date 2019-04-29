// ECE556 - Copyright 2014 University of Wisconsin-Madison.  All Rights Reserved.

#include "ece556.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "limits.h"
#include "float.h"
#include <stdlib.h>
#include <vector>
#include <set>
#include <stack>
#include <queue>
#include <iostream>

using namespace std;

typedef unsigned long long ulonglong;

typedef struct {
  point p1;
  point p2;
  ulonglong weight;
} edge;

typedef struct {
  bool operator()(const edge& e1, const edge& e2) {
    return e1.weight > e2.weight;
  }
} less_PQ;

typedef struct {
  bool operator()(const point& p1, const point& p2) {
    return (p1.x < p2.x) || ((p1.x == p2.x) && (p1.y < p2.y));
  }
} less_point;

int panic(const char* mssg) {
  printf("panic: %s\n", mssg);
  return 0;
}

// Gets the distance between two points
double getDistance(const point *p1, const point *p2){
  
  double dy = abs( (double) ((p1->y) - (p2->y)) );  // Vertical separation between points
  double dx = abs( (double) ((p1->x) - (p2->x)) );  // Horizontal separation between points
  return sqrt(dy*dy + dx*dx);          // Pythagorean theorem to find distance between points
}

/*int netDecompDjikstra(routingInst *rst) {
  // Assumes net_decomp_order has already been allocated.
  for(int i = 0; i < rst->numNets; i++) {// Loops for each net.
    // Reorder each net.
    pair<int, point> start = {0, {rst->nets[0].x rst->nets[0].y}};
    for(int j = 0; j < rst->nets[i].numPins; j++) {// Loops for each pin in given net.
    }
  }
}*/

// Reorders the pins within each net
/*int netDecomp(routingInst *rst) {
  // For each net
  for (int n = 0; n < rst->numNets; n++) {
    
  std::vector<point> unRouted;
  std::vector<point> routed;
    
  // Add all pins to the "unRouted" list
  for (int p = 0; p < rst->nets[n].numPins; p++) {
    unRouted.push_back(rst->nets[n].pins[p]);
  }
    
  // Move bottom pin out of unRouted list and into routed list
  int bottom = 0;
  int miny = INT_MAX;
  for (unsigned int unro = 0; unro < unRouted.size(); unro++) {
    if (unRouted[unro].y < miny ) {
      bottom = unro;
      miny = unRouted[unro].y;
    }
  }
  routed.push_back( unRouted[bottom] );
  unRouted.erase( unRouted.begin() + bottom);
    
  // Store the closest neighbor
  int bestPinIndex;
  double bestPinDist;
    
  // While there are still unRouted pins
  while (!unRouted.empty()) {
    bestPinIndex = 0;
    bestPinDist = DBL_MAX;
        
    // Routing from this routed pin
    point source = routed.back();
        
    // Find the best unRouted pin to route
    for (unsigned int unro = 0; unro < unRouted.size(); unro++) {
      // is this the new best?
      double thisDist = getDistance(&source, &(unRouted[unro]));
      if ( thisDist < bestPinDist) {
	bestPinIndex = unro;
	bestPinDist = thisDist;
      }
    }
        
    // Move the best pin into routed list
    routed.push_back(unRouted[bestPinIndex]);
        
    // Remove the best pin from the unRouted list
    unRouted.erase(unRouted.begin() + bestPinIndex);
  }
    
    // Fill in the pins for this net from the routed list
    for (int p = 0; p < rst->nets[n].numPins; p++) {
      rst->nets[n].pins[p] = routed[p];
    }
  }
  
  return 1;
}*/



/*
  Finds the edge number that is between two adjacent points.
  This will panic if the two points aren't adjacent.
*/
int findEdge(const point *p1, const point *p2, const int gx, const int gy){
  
  if(p1 == NULL){
    panic("findEdge: p1 pointer was NULL");
    return -1;
  }

  if(p2 == NULL){
    panic("findEdge: p2 pointer was NULL");
    return -1;
  }

  int distance = (int) getDistance(p1, p2);

  if(distance != 1){
    panic("findEdge: Not an edge -- distance between points");
    return -1;
  }

  // Check if this edge is Vertical or Horizontal
  if( (p1->x) == (p2->x) ){
    // This is a vertical edge
    int x = p1->x;                    // x coordinate of edge
    int y = std::min(p1->y,p2->y);    // lower y coordinate of edge

    return ((gx - 1) * gy) + (gx * y) + x;  // return the edge number
  } else {
    // This is a horizontal edge
    int y = p1->y;                    // y coordinate of edge
    int x = std::min(p1->x,p2->x);    // left-most x coordinate of edge

    return ((gx - 1) * y) + x;        // return the edge number
  }

}

int readBenchmark(const char *fileName, routingInst *rst){
  /*********** TO BE FILLED BY YOU **********/  
  if(fileName == NULL || rst == NULL)
    return panic("readBenchmark: null argument");
  // Used to read file.
  FILE* fp;
  char* ln_str = NULL;
  size_t ln_sz = 0;
  // Used for processing file.
  char* ln_cpy = NULL;
  char* word = NULL;
  // Variables for routingInst
  int numNets;
  int numEdges;
  int cap;

  fp = fopen(fileName, "r");
  if(fp == NULL)
    return panic("file could not be opened.");

  // Retreive grid dimensions.
  if(getline(&ln_str, &ln_sz, fp) == -1)
    panic("grid: couldn't read input file.");
  ln_cpy = (char*)malloc(sizeof(char) * ln_sz);
  ln_cpy = strcpy(ln_cpy, ln_str);
  word = strtok(ln_cpy, "\t");
  if(strcmp(word, "grid") != 0)
    return panic("Couldn't find word grid.");
  word = strtok(NULL, " ");
  rst->gx = atoi(word);
  word = strtok(NULL, "\n");
  rst->gy = atoi(word);
  if(strtok(NULL, " ") != NULL)
    return panic("First line was too long.");
  free(ln_cpy);
  ln_cpy = NULL;

  // Retreive capacity.
  if(getline(&ln_str, &ln_sz, fp) == -1)
    return panic("capacity: couldn't read input file.");
  ln_cpy = (char*)malloc(sizeof(char) * ln_sz);
  ln_cpy = strcpy(ln_cpy, ln_str);
  word = strtok(ln_cpy, "\t");
  if(strcmp(word, "capacity"))
    return panic("Couldn't find word capacity");
  word = strtok(NULL, "\n");
  cap = atoi(word);// Extra variable reduces runtime when it's used later.
  rst->cap = cap;
  if(strtok(NULL, "\t") != NULL)
    return panic("First line was too long.");
  free(ln_cpy);
  ln_cpy = NULL;

  // Retreive number of nets.
  if(getline(&ln_str, &ln_sz, fp) == -1)
    return panic("num nets: couldn't read input file.");
  ln_cpy = (char*)malloc(sizeof(char) * ln_sz);
  ln_cpy = strcpy(ln_cpy, ln_str);
  word = strtok(ln_cpy, " ");
  if(strcmp(word, "num"))
    return panic("Couldn't find word num.");
  word = strtok(NULL, " ");
  if(strcmp(word, "net"))
    return panic("Couldn't find word net.");
  word = strtok(NULL, "\n");
  numNets = atoi(word);
  rst->numNets = numNets;// Extra variable to keep info on stack for quick access.
  if(strtok(NULL, " ") != NULL)
    return panic("First line was too long.");
  free(ln_cpy);
  ln_cpy = NULL;

  // Allocate arrays.
  rst->nets = (net*)malloc(sizeof(net) * numNets);
  numEdges = (rst->gx - 1) * rst->gy + (rst->gy - 1) * rst->gx;
  rst->numEdges = numEdges;
  rst->edgeCaps = (int*)malloc(sizeof(int) * numEdges);
  rst->edgeUtils = (int*)malloc(sizeof(int) * numEdges);
  for(int i = 0; i < numEdges; i++) {// Initialize edgeCaps to default value.
    rst->edgeCaps[i] = cap;
  }
  for(int i = 0; i < numEdges; i++) {// Initialize edgeUtils to 0.
    rst->edgeUtils[i] = 0;
  }

  for(int i = 0; i < numNets; i++) {
    if(getline(&ln_str, &ln_sz, fp) == -1)
      return panic("num pins: couldn't read input file.");
    ln_cpy = (char*)malloc(sizeof(char) * ln_sz);
    ln_cpy = strcpy(ln_cpy, ln_str);
    word = strtok(ln_cpy, " ");
    const int net_num = atoi(&word[1]);
    rst->nets[i].id = net_num;
    word = strtok(NULL, "\n");
    const int num_pins = atoi(word);
    rst->nets[i].numPins = num_pins;
    point* pins = (point*)malloc(sizeof(point) * num_pins);
    if(strtok(NULL, " ") != NULL)
      return panic("num pins: line is too long.");
    free(ln_cpy);
    ln_cpy = NULL;

    for(int j = 0; j < num_pins; j++) {
      if(getline(&ln_str, &ln_sz, fp) == -1)
        return panic("pin traversal: couldn't read input file.");
      ln_cpy = (char*)malloc(sizeof(char) * ln_sz);
      ln_cpy = strcpy(ln_cpy, ln_str);
      word = strtok(ln_cpy, "\t");
      pins[j].x = atoi(word);
      word = strtok(NULL, "\n");
      pins[j].y = atoi(word);
      if(strtok(NULL, "\t") != NULL)
        return panic("pin traversal: line is too long.");
      free(ln_cpy);
      ln_cpy = NULL;
    }
    rst->nets[i].pins = pins;
  }

  // Read blockages.
  if(getline(&ln_str, &ln_sz, fp) == -1)
    return panic("blockages: couldn't read lines.");
  const int num_blkgs = atoi(ln_str);
  for(int i = 0; i < num_blkgs; i++) {
    if(getline(&ln_str, &ln_sz, fp) == -1)
      return panic("blockages: couldn't read block points.");
    ln_cpy = (char*)malloc(sizeof(char) * ln_sz);
    ln_cpy = strcpy(ln_cpy, ln_str);
    point* p1 = (point*)malloc(sizeof(point));
    point* p2 = (point*)malloc(sizeof(point));
    word = strtok(ln_cpy, " ");
    p1->x = atoi(word);
    word = strtok(NULL, "\t");
    p1->y = atoi(word);
    word = strtok(NULL, " ");
    p2->x = atoi(word);
    word = strtok(NULL, " ");
    p2->y = atoi(word);
    word = strtok(NULL, "\n");
    const int new_cap = atoi(word);
    const int edge = findEdge(p1, p2, rst->gx, rst->gy);
    if(edge < 0)
      return panic("blockages: invalid edge.");
    rst->edgeCaps[edge] = new_cap;
    if(strtok(NULL, " ") != NULL)
      return panic("blockages: line is too long");
    free(ln_cpy);
  }
  fclose(fp);
  return 1;
}


/*
Routes a segment (a straight horiz or vert line) that connects the two given points.
Increments the edge utilizations as it goes.
The two points MUST SHARE A COMMON X OR Y COORDINATE! (otherwise panic and exit...)
*/
int routeSeg(point p1, point p2, int netNum, routingInst *rst) {
  
  // Check that the two points form a horiz or vert line
  if ( (p1.x != p2.x) && (p1.y != p2.y) ) {
    // These two points should not have been given
    panic("routeSeg: The two points you gave me were not a horiz or vert pair!");
    exit(EXIT_FAILURE);
  }

  // segment that will be pushed to the vector
  segment seg;

  // define seg's points
  seg.p1 = p1;
  seg.p2 = p2;

  // will the segment be vertical or horizontal?
  if (p1.y == p2.y) {
    // Horizontal Pair

    // Make p1 the left point
    if (p1.x > p2.x) {
      point temp = p1;
      p1 = p2;
      p2 = temp;
    }

    // Separation between the two points
    int dx = p2.x - p1.x;

    // Now we know how many edges there are
    seg.numEdges = dx;

    // Allocate heap memory for the edges array
    seg.edges = (int*)malloc(sizeof(int) * seg.numEdges);
    if (seg.edges == NULL) {
      panic("routeSeg: Failure allocating heap memory (horizontal)");
      exit(EXIT_FAILURE);
    }

    point left, right;
    
    // Fill in the edges and increment edgeUtils
    for (int e = 0; e < seg.numEdges; e++) {
      left.x = p1.x + e;
      left.y = p1.y;

      right.x = left.x + 1;
      right.y = p1.y;
      
      int edge = findEdge( &left, &right, rst->gx, rst->gy);
      seg.edges[e] = edge;
      rst->edgeUtils[edge]++;
    }
  } else {
    // Vertical Pair

    // Make p1 the bottom point
    if (p1.y > p2.y){
      point temp = p1;
      p1 = p2;
      p2 = temp;
    }

    // Separation between the two points
    int dy = p2.y - p1.y;

    // Now we know how many edges there are
    seg.numEdges = dy;

    // Allocate heap memory for the edges array
    seg.edges = (int*)malloc(sizeof(int) * seg.numEdges);
    if (seg.edges == NULL) {
      panic("routeSeg: Failure allocating heap memory (vertical)");
      exit(EXIT_FAILURE);
    }
    
    point bottom, top;
    
    // Fill in the edges and increment edgeUtils
    for (int e = 0; e < seg.numEdges; e++) {
      bottom.x = p1.x;
      bottom.y = p1.y + e;

      top.x = p1.x;
      top.y = bottom.y + 1;
      
      int edge = findEdge( &bottom, &top, rst->gx, rst->gy);
      seg.edges[e] = edge;
      rst->edgeUtils[edge]++;
    }
  }  
  // Push the segment onto the vector
  rst->nets[netNum].nroute.segments.push_back(seg);

  return 1;
}

int routePinsBasic(point p1, point p2, int netNum, routingInst *rst){
  // Check if this is a diagonal pair
  if ( (p1.x != p2.x) && (p1.y != p2.y) ) {
    // It's a diagonal pair
    point corner;
    corner.x = p2.x;
    corner.y = p1.y;
    routeSeg(p1, corner, netNum, rst);
    routeSeg(corner, p2, netNum, rst);
  } else {
    // They could be horizontal or vertical
    routeSeg(p1, p2, netNum, rst);
  }
	
  return 1;
}

/*
Computes the best route between the two given points in the given net.
*/
int routePins(point p1, point p2, int netNum, routingInst *rst) {
  
  if((p1.x == p2.x) && (p1.y == p2.y))
    panic("Same point for p1 and p2.");
  set<point, less_point> P_set;
  stack<edge> RP;// Used for retrace path.
  const int width = rst->gx;
  const int height = rst->gy;
  const int capacity = rst->cap;
  const int coeff = (capacity / 8) + 1;// Will be multiplied by length to destination for heuristic.
  point new_point;
  edge new_edge;
  ulonglong Q2_weight = 0;// Edge weight.
  int edge_num;// Edge number.
  priority_queue<edge, vector<edge>, less_PQ> Q2;// Queue that will be hold points for Q2.
  pair<int, point> start_point = {0, p1};// Point that edges will be taken from for Q2.
  
  while(true) {
    // Add all neighboring edges of p1 to Q2.
    // Add right edge.
    if((start_point.second.x >= 0) && (start_point.second.x + 1 < width)) {
      new_point = {start_point.second.x + 1, start_point.second.y};
      if(P_set.find(new_point) == P_set.end()) {
        edge_num = findEdge(&start_point.second, &new_point, width, height);
        Q2_weight = start_point.first + capacity - rst->edgeCaps[edge_num] + rst->edgeUtils[edge_num];// Heuristic for weight.
        if(Q2.size() == 0)
          Q2_weight += (abs(p2.x - new_point.x) + abs(p2.y - new_point.y)) * coeff;// Heuristic for estimated weight to destination.
	else if(p2.x - start_point.second.x > 0)
          Q2_weight -= coeff;
	else if(p2.x - start_point.second.x < 0)
          Q2_weight += coeff;
        new_edge = {start_point.second, new_point, Q2_weight};// pair to add.
        Q2.push(new_edge);
        //P_set.insert(new_point);//FIXME
      }
    }
    // Add top edge.
    if((start_point.second.y >= 0) && (start_point.second.y + 1 < height)) {
      new_point = {start_point.second.x, start_point.second.y + 1};
      if(P_set.find(new_point) == P_set.end()) {
        edge_num = findEdge(&start_point.second, &new_point, width, height);
        Q2_weight = start_point.first + capacity - rst->edgeCaps[edge_num] + rst->edgeUtils[edge_num];// Heuristic for weight.
        if(Q2.size() == 0)
          Q2_weight += (abs(p2.x - new_point.x) + abs(p2.y - new_point.y)) * coeff;// Heuristic for estimated weight to destination.
	else if(p2.y - start_point.second.y > 0)
          Q2_weight -= coeff;
	else if(p2.y - start_point.second.y < 0)
          Q2_weight += coeff;
        new_edge = {start_point.second, new_point, Q2_weight};// pair to add.
        Q2.push(new_edge);
        //P_set.insert(new_point);
      }
    }
    // Add left edge.
    if((start_point.second.x > 0) && (start_point.second.x < width)) {
      new_point = {start_point.second.x - 1, start_point.second.y};
      if(P_set.find(new_point) == P_set.end()) {
        edge_num = findEdge(&start_point.second, &new_point, width, height);
        Q2_weight = start_point.first + capacity - rst->edgeCaps[edge_num] + rst->edgeUtils[edge_num];// Heuristic for weight.
        if(Q2.size() == 0)
          Q2_weight += (abs(p2.x - new_point.x) + abs(p2.y - new_point.y)) * coeff;// Heuristic for estimated weight to destination.
	else if(p2.x - start_point.second.x > 0)
          Q2_weight += coeff;
	else if(p2.x - start_point.second.x < 0)
          Q2_weight -= coeff;
        new_edge = {start_point.second, new_point, Q2_weight};// pair to add.
        Q2.push(new_edge);
        //P_set.insert(new_point);
      }
    }
    // Add bottom edge.
    if((start_point.second.y > 0) && (start_point.second.y < height)) {
      new_point = {start_point.second.x, start_point.second.y - 1};
      if(P_set.find(new_point) == P_set.end()) {
        edge_num = findEdge(&start_point.second, &new_point, width, height);
        Q2_weight = start_point.first + capacity - rst->edgeCaps[edge_num] + rst->edgeUtils[edge_num];// Heuristic for weight.
        if(Q2.size() == 0)
          Q2_weight += (abs(p2.x - new_point.x) + abs(p2.y - new_point.y)) * coeff;// Heuristic for estimated weight to destination.
	else if(p2.x - start_point.second.x > 0)
          Q2_weight += coeff;
	else if(p2.x - start_point.second.x < 0)
          Q2_weight -= coeff;
        new_edge = {start_point.second, new_point, Q2_weight};// pair to add.
        Q2.push(new_edge);
        //P_set.insert(new_point);
      }
    }
    // Add point with new Q2 edges to Q3.
    //Q3.push({Q2.top().p1, Q2.top().p2, Q2.top().weight});
    while(P_set.find(Q2.top().p2) != P_set.end())
      Q2.pop();
    RP.push(Q2.top());
    P_set.insert(Q2.top().p2);
    //printf("Added edge to Q3.\n");//FIXME
    // If it is the last point, use retrace path to route all remaining points.
    if((Q2.top().p2.x == p2.x) && (Q2.top().p2.y == p2.y)) {
      //printf("Found one point.\n");//FIXME
      break;
    }
    // Find lowest cost edge from Q2.
    // Add edge to retrace path, and set point as start to add to Q3 if it's not p2.
    start_point = {Q2.top().weight, Q2.top().p2};
    Q2.pop();
  }

  //TODO: Figure out if point should be added after it's in Q2 or Q3. Could reduce runtime.
  if((RP.top().p1.x == p1.x) && (RP.top().p1.y == p1.y)) {
    routeSeg(RP.top().p1, RP.top().p2, netNum, rst);
    return 1;
  }

  bool done = false;
  point end_point = RP.top().p2;
  while(!done) {
    bool horizontal = true;
    bool vertical = true;
    new_edge = RP.top();
    RP.pop();
    while((RP.top().p2.x != new_edge.p1.x) || (RP.top().p2.y != new_edge.p1.y)) {
      if(RP.size() == 0)
        panic("About to pop empty stack!");
      RP.pop();
    }
    //printf("Found one edge to route.\n");//FIXME
    horizontal = (RP.top().p1.x == end_point.x);
    vertical = (RP.top().p1.y == end_point.y);
    if((!horizontal && !vertical) || ((RP.top().p1.x == p1.x) && (RP.top().p1.y == p1.y))) {
      if((RP.top().p1.x == p1.x) && (RP.top().p1.y == p1.y)) {
        if(!horizontal && !vertical) {
          routeSeg(RP.top().p2, end_point, netNum, rst);
          routeSeg(RP.top().p1, RP.top().p2, netNum, rst);
	}
	else {
          routeSeg(RP.top().p1, end_point, netNum, rst);
	}
        done = true;
      }
      else {
        //printf("About to call routeSeg.\n");//FIXME
        routeSeg(RP.top().p2, end_point, netNum, rst);
        //printf("Called routeSeg.\n");//FIXME
        end_point = RP.top().p2;
      }
    }
  }
	

  // Finished
  return 1;
}

// Returns the total overflow in the routing instance
int getOverflow(routingInst *rst){
  // Record overflow
  int total_overflow = 0;
  int overflow_here;
  
  // For each edge on the grid
  for (int e = 0; e < rst->numEdges; e++){
    // Get the overflow at this edge
    overflow_here = rst->edgeUtils[e] - rst->cap;
    
    // Check if this edge has overflow or not
    if(overflow_here > 0){
      // add this overflow to the total
      total_overflow = total_overflow + overflow_here;
    }
  }

  // Return the total overflow
  return total_overflow;
}

// Will route pins by making a minimum spanning tree.
void routeNetDcmp(routingInst *rst, int netNum, int ripUp) {
  bool tree_done = false;
  const uint num_points = rst->nets[netNum].numPins;
  pair<int, point> start = {0, rst->nets[netNum].pins[0]};
  priority_queue<edge, vector<edge>, less_PQ> PQ;
  vector<edge> tree;
  set<point, less_point> P_set;

  P_set.insert(start.second);

  while(!tree_done) {
    // Add edges for every point to connect to start.
    for(int p = 0; p < rst->nets[netNum].numPins; p++) {
      const point new_p = rst->nets[netNum].pins[p];
      if(P_set.find(new_p) == P_set.end()) {
        ulonglong weight = abs(start.second.x - new_p.x) + abs(start.second.y - new_p.y);
        PQ.push({start.second, new_p, weight});
      }
    }
    while(P_set.find(PQ.top().p2) != P_set.end())// Pop all points already in minimum spanning tree.
      PQ.pop();
    tree.push_back(PQ.top());
    P_set.insert(PQ.top().p2);
    start = {PQ.top().weight, PQ.top().p2};
    PQ.pop();
    if(P_set.size() == num_points) {
      tree_done = true;
      break;
    }
  }
  for(uint e = 0; e < tree.size(); e++) {
    if(ripUp)
      routePins(tree[e].p1, tree[e].p2, netNum, rst);
    else
      routePinsBasic(tree[e].p1, tree[e].p2, netNum, rst);
  }
}

int solveRoutingBasic(routingInst *rst, int netDcmp, int ripUp){
 
  if(netDcmp) {
    for(int n = 0; n < rst->numNets; n++) {
      routeNetDcmp(rst, n, ripUp);
    }
  }
  else {  
    // For each net
    for (int n = 0; n < rst->numNets; n++) {
      // For all but the last pin in this net
      for (int p = 0; p < ((rst->nets[n].numPins) - 1); p++) {
        // Route the pair (thisPin, nextPin)
        routePinsBasic((rst->nets[n].pins[p]), (rst->nets[n].pins[p + 1]), n, rst);
      }
    }
  }

  return 1;
}

/*
int solveRouting(routingInst *rst, int netDcmp){//FIXME: Possibly remove this function?
  
  if(netDcmp) {
    for(int n = 0; n < rst->numNets; n++) {
      routeNetDcmp(rst, n);
    }
  }
  else {
    // For each net
    for (int n = 0; n < rst->numNets; n++) {
      // For all but the last pin in this net
      for (int p = 0; p < ((rst->nets[n].numPins) - 1); p++) {
        // Route the pair (thisPin, nextPin)
        routePins((rst->nets[n].pins[p]), (rst->nets[n].pins[p + 1]), n, rst);
      }
    }
  }
  return 1;
}*/

// This uses a shuffled order
int solveRoutingRand(int *order, routingInst *rst, int netDcmp){
  
  if(netDcmp) {
    for(int n = 0; n < rst->numNets; n++) {
      routeNetDcmp(rst, n, 1);
    }
  }
  else {
    // For each net
    for (int n = 0; n < rst->numNets; n++) {
      // For all but the last pin in this net
      for (int p = 0; p < ((rst->nets[order[n]].numPins) - 1); p++) { 
        // Route the pair (thisPin, nextPin)
        routePins((rst->nets[order[n]].pins[p]), (rst->nets[order[n]].pins[p + 1]), order[n], rst);
      }
    }
  }
  return 1;
}

int writeOutput(const char *outRouteFile, routingInst *rst){
  /*********** TO BE FILLED BY YOU **********/
  FILE* fp = NULL;
  if((fp = fopen(outRouteFile, "w")) == NULL)
    panic("Could not open file.");
  const int numNets = rst->numNets;
  for(int n = 0; n < numNets; n++) {// Iterates for each net.
    fprintf(fp, "n%d\n", n);
    for(const segment& seg : rst->nets[n].nroute.segments) {// Iterates for each segment or edge.
      fprintf(fp, "(%d,%d)-(%d,%d)\n", seg.p1.x, seg.p1.y, seg.p2.x, seg.p2.y);
    }
    fprintf(fp, "!\n");
  }
  fclose(fp);
  return 1;
}

/*
int releaseNetsOnly(routingInst *rst){
	const int numNets = rst->numNets;
	for(int n = 0; n < numNets; n++) {
		free(rst->nets[n].pins);
		for(const segment& seg : rst->nets[n].nroute.segments) {
			free(seg.edges);
		}
	}
	return 1;
}
*/


/*
	Shuffles the order of the nets in the routing instance
*/
int *makeOrderArray(routingInst *rst){
	// Make the array
	int *order = (int *)malloc(sizeof(int) * (rst->numNets));
	for (int i = 0; i < rst->numNets; i++) {
		order[i] = i;
	}
	/*
	// Rotate the array left
	for (int i = 0; i < (rst->numNets - 1); i++) {
		order[i] = order[i + 1];
	}*/
	
	//order[rst->numNets - 1] = 0;
	
	
	// Shuffle the array
	int n = rst->numNets;
	
	if (n > 1) 
    {
		
        int i;
        for (i = 0; i < n - 1; i++) 
        {
          int j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = order[j];
          order[j] = order[i];
          order[i] = t;
        }
    }
	

	return order;
}
/*
int unShuffleNets(routingInst *rst){
	
	net *Array = rst->nets;
	net temp;
	int i, j;
	int Size = rst->numNets;
	
	for (i = 0; i < Size; i++)
	{
		for (j = i + 1; j < Size; j++)
		{
			if(Array[i].id > Array[j].id)
			{
				temp = Array[i];
				Array[i] = Array[j];
				Array[j] = temp;
			}
			
		}
	}
	
	return 1;
}
*/
int release(routingInst *rst){
  /*********** TO BE FILLED BY YOU **********/
  const int numNets = rst->numNets;
  for(int n = 0; n < numNets; n++) {
    free(rst->nets[n].pins);
    for(const segment& seg : rst->nets[n].nroute.segments) {
      free(seg.edges);
    }
    // TODO: Add rst->nets[n].nroute.segments.clear() here if necessary
  }
  free(rst->edgeCaps);
  free(rst->edgeUtils);
  return 1;
} 
