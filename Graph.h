#ifndef GRAPH_H_
#define GRAPH_H_
#include "Motif.h"
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
using namespace std;

class Graph{
public:
	int NodeSize;
	int EdgeSize;
	int LayerSize;
	//aggregate graph
    vector<dynamic_bitset<> > adjacencyList;
	int **adjacencyMatrix;
	int k;
	bool cannotSwap;
	vector<dynamic_bitset<> >edgeLayers; //for each edge, store the set of layers that have this edge
	vector<pair<int, int>> edgeIntSet; //for each edge, record the mapping id
	vector<pair<string, string>> edgeset; //for each edge, record original node name
	
public:
	Graph(int layerNum, int nodeNum, string networkName, double thre, int swapNumber, int swapMethod): LayerSize(layerNum), NodeSize(nodeNum), adjacencyList(nodeNum, dynamic_bitset<>(nodeNum)){
		cannotSwap = false;
		k = ceil(layerNum*thre);
		if(swapMethod == 1 && swapNumber > 0) networkName = randomShuffleEdge(networkName, swapNumber);
		else if(swapMethod == 2 && swapNumber > 0){
			networkName = nonAggregateShuffleEdge(networkName, swapNumber);
			if(networkName == "") cannotSwap = true;
		}else if(swapMethod == 3 && swapNumber > 0){
			networkName = aggregateShuffleEdge(networkName, swapNumber);
		}else if(swapMethod == 4 && swapNumber > 0){
			networkName = aggNonAggregateShuffleEdge(networkName, swapNumber);
		}

		if(!cannotSwap){
			adjacencyMatrix = new int*[nodeNum];
			for(int i=0;i<nodeNum;i++){
				adjacencyMatrix[i] = new int[nodeNum]();
			}
			readNetworkFile(networkName);
		}
	}

    ~Graph(){
		for(int i=0;i<NodeSize;i++) {
			 delete[] adjacencyMatrix[i];
		}
		delete[] adjacencyMatrix;
		adjacencyList.clear();
	}

    void findMotifs(vector<Motif*>& pattern, int caseN);
	string randomShuffleEdge(string networkName, int swapNumber);
	string nonAggregateShuffleEdge(string networkName, int swampNumber);
	string aggregateShuffleEdge(string networkName, int swapNumber);
	string aggNonAggregateShuffleEdge(string networkName, int swapNumber);
	string motifShuffle(vector<Motif*>& motifs, string networkName, int swapNumber);
private:
    void readNetworkFile(string networkName);
};

#endif
