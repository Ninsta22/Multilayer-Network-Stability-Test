#include "Graph.h"
#include <iostream>
#include <chrono>
#include <map>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <cstdlib>
Graph* g;

int F2count(vector<Motif*>& pattern, vector<int>& selectEdge, dynamic_bitset<>& networks){
	int countF2 = 0;
	int size = pattern.size();
	if (size == 0) {
		return countF2;
	}

	//build overlap graph
	vector<dynamic_bitset<> > madjacencyList(size, dynamic_bitset<>(size));
	vector<dynamic_bitset<> > layerAdj(g->LayerSize, dynamic_bitset<>(size));
	dynamic_bitset<>::size_type it, itt;

	for(int i=0;i<size;i++){
		Motif* a=pattern[i];

		for (int j = 0; j < g->LayerSize; j++) {
			if (a->layers[j] == 1) {
				layerAdj[j][i] = 1;
			}
		}

		for(int j=i+1;j<size;j++){
			Motif* b = pattern[j];
			dynamic_bitset<> tmp = b->edges;
			tmp &= a->edges;
			if (tmp.any()) {
				//			cout<<"neighbors: "<<j<<endl;
				madjacencyList[i][j] = 1;
				madjacencyList[j][i] = 1;
			}
		}
	}

	dynamic_bitset<> neighbors(size);
	dynamic_bitset<> embeddings(size);
	embeddings.set();
	networks.set();

	while (embeddings.any()) {
		double curMinLoss = embeddings.count()+1;
        int index = -1;
		//calculate embeddings's loss
		it=embeddings.find_first();
		while(it!=dynamic_bitset<>::npos){
			Motif* b=pattern[it];
			dynamic_bitset<> curNetworks = b->layers;

			neighbors = madjacencyList[it];
			double loss = neighbors.count();
			neighbors[it] = 1;
			dynamic_bitset<> rmLayers = curNetworks;
			rmLayers ^= networks;

			if (rmLayers.any() && loss <= curMinLoss) {
				dynamic_bitset<> others(size);
				itt = rmLayers.find_first();
				while(itt != dynamic_bitset<>::npos) {
					others |= layerAdj[itt];
					itt = rmLayers.find_next(itt);
				}

				dynamic_bitset<> tmp = neighbors;
				tmp &= others;
				others ^= tmp;

				itt = others.find_first();
				while(itt!=dynamic_bitset<>::npos && loss <= curMinLoss){
					Motif* c = pattern[itt];
					dynamic_bitset<> tmp = c->layers;
					tmp &= curNetworks;
					if (tmp.count() < g->k) {
						neighbors[itt] = 1;
						loss += 1;
					}
					else {
						loss += 0.1*(1 - tmp.count()*1.0/c->layers.count());
					}
					itt=others.find_next(itt);
				}
			}

			if (loss < curMinLoss) {
				curMinLoss = loss;
				index = it;
			}
			b->removes = neighbors;
			it=embeddings.find_next(it);
		}
        
		//select the candidates with least loss
		Motif* a = pattern[index];

		//		cout<<"select "<<a->id<<" with loss "<<a->degree<<endl;
		countF2++;
		networks &= a->layers;
		for (int i = 0; i < g->EdgeSize; i++) {
			if (a->edges[i] == 1) {
				selectEdge.push_back(i);
			}
		}

		//delete a and its neighbors
		neighbors = a->removes;
		for (int i = 0; i < g->LayerSize; i++) {
			dynamic_bitset<> tmp = layerAdj[i];
			tmp &= neighbors;
			layerAdj[i] ^= tmp;
		}

		it = embeddings.find_first();
		while(it!=dynamic_bitset<>::npos){
			Motif* b = pattern[it];
			if (neighbors[it] == 1) {
				embeddings[it] = 0;
				delete pattern[it];
			}
			else {
				b->layers &= networks;
				dynamic_bitset<> tmp = madjacencyList[it];
				tmp &= neighbors;
				madjacencyList[it] ^= tmp;
			}
			it=embeddings.find_next(it);
		}
	}

	neighbors.clear();
	pattern.clear();
	madjacencyList.clear();
	return countF2;
}

int main(int argc, char *argv[]){
	//parameters: which motif, numberOfLayers, numberOfNodes, folderName, threshold
	int type = stoi(argv[1]);
	int numberLayers = stoi(argv[2]);
	int numberNodes = stoi(argv[3]);
	string networkName = argv[4];
	double threshold = stod(argv[5]);

	auto f2_start = chrono::steady_clock::now();
	g = new Graph(numberLayers, numberNodes, networkName, threshold);
	//find all possible embeddings
	vector<Motif*> pattern;
	g->findMotifs(pattern,type);
    
    //find F2 count, selectEdges stores the edges of final motif set, and networks stores the layers which have final motif set
	vector<int> selectEdge;
	dynamic_bitset<> networks(g->LayerSize);
	int numEmbeddings = F2count(pattern, selectEdge, networks);
	auto f2_end = chrono::steady_clock::now();
	cout<<"Number of embedding\tRunning time(ms)"<<endl;
	cout<<numEmbeddings<<"\t"<<chrono::duration_cast<chrono::milliseconds>(f2_end - f2_start).count()<<endl;
	//print layers of motifs
    cout<<"The motifs are on layers (starting from layer 1):"<<endl;
    for (int i = 0; i < g->LayerSize; i++) {
    	if (networks[i] == 1) {
    		cout<<"Layer "<<i + 1<<endl;
    	}
    }

    //save edges of final motif set
    ofstream finalRes;
    finalRes.open("motif_edges.txt");
    for (int i = 0; i < selectEdge.size(); i++) {
        finalRes<<g->edgeset[selectEdge[i]].first<<"\t"<<g->edgeset[selectEdge[i]].second<<endl;
    }
    finalRes.close();
	delete g;

	
}
