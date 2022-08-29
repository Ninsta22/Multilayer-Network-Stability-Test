#ifndef GRAPH_H_
#define GRAPH_H_
#include "Motif.h"
#include <vector>
#include <cmath>
#include <string>
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
	vector<dynamic_bitset<> >edgeLayers; //for each edge, store the set of layers that have this edge
	vector<pair<int, int>> edgeIntSet; //for each edge, record the mapping id
	vector<pair<string, string>> edgeset; //for each edge, record original node name

public:
	Graph(int layerNum, int nodeNum, string networkName, double thre): LayerSize(layerNum), NodeSize(nodeNum), adjacencyList(nodeNum, dynamic_bitset<>(nodeNum)){
        k = ceil(layerNum*thre);
		adjacencyMatrix = new int*[nodeNum];
		for(int i=0;i<nodeNum;i++)
		    adjacencyMatrix[i] = new int[nodeNum]();
		readNetworkFile(networkName);
	}

    ~Graph(){
		for(int i=0;i<NodeSize;i++) {
			 delete[] adjacencyMatrix[i];
		}
		delete[] adjacencyMatrix;
		adjacencyList.clear();
	}

    void findMotifs(vector<Motif*>& pattern, int caseN);
private:
    void readNetworkFile(string networkName);
};

#endif
