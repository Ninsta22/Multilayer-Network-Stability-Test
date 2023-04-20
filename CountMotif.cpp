#include "Graph.h"
#include <iostream>
#include <chrono>
#include <map>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <cstdlib>
Graph* g;


vector<int> getNodeLayerInfo(string networkName){
	//Count number of layers and nodes in file
	ifstream in_stream;
	in_stream.open(networkName);
	int layerCounter = 0;
	int largestNode = -1;
	string str;
	vector<int> returnInformation;

	while(!in_stream.eof()){
		getline(in_stream, str);
		if(str.size()>0){
			if(str[0] == '-') layerCounter++;
			else{
				bool flag = true;
				string firstNumber = "";
				string secondNumber = "";

				for(int i = 0; i < str.length(); i++){
					if(str[i] == '\t') flag = false;
					else if(flag) firstNumber = firstNumber + str[i];
					else secondNumber = secondNumber + str[i];
				}

				int first = stoi(firstNumber);
				int second = stoi(secondNumber);

				if(largestNode < first) largestNode = first;
				if(largestNode < second) largestNode = second;
			}
		}
	}

	in_stream.close();
	returnInformation.push_back(layerCounter);
	returnInformation.push_back(largestNode+1);
	return returnInformation;
}

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
	//parameters: which motif, folderName, swapping method
	
	/*
	int type = stoi(argv[1]);
	string networkName = argv[2];
	double threshold = stod(argv[3]);
	int swapNumber = stoi(argv[4]);*/



	int type = stoi(argv[1]);
	string startingNetworkName = argv[2];
	int swapMethod = stoi(argv[3]);
	string csvName = "";
	string motifName = "";
	
	if(swapMethod == 1) csvName = "pure_random";
	else if(swapMethod == 2) csvName = "non_aggregate";
	else if(swapMethod == 3) csvName = "aggregate";
	else if(swapMethod == 4) csvName = "non_aggregate_with_aggregate";
	else if(swapMethod == 5) csvName = "stochastic";
	else if(swapMethod == 6) csvName = "total";

	ofstream totalCSV;
    totalCSV.open(csvName+"_results.csv");
	totalCSV << "swapType, motifType, threshold, totalSwaps, Number of Motifs, Time(ms)"<<endl;

	vector<string> swapFileNames;
	swapFileNames.push_back(startingNetworkName); swapFileNames.push_back("10swapNetwork.txt"); swapFileNames.push_back("15swapNetwork.txt");
	swapFileNames.push_back("25swapNetwork.txt"); swapFileNames.push_back("50swapNetwork.txt");
	vector<double> testingNumbers;
	//testingNumbers.push_back(0.2); testingNumbers.push_back(0.4); testingNumbers.push_back(0.6);
	testingNumbers.push_back(0.4);
	vector<int> listOfSwapNumbers;
	listOfSwapNumbers.push_back(0); listOfSwapNumbers.push_back(10); listOfSwapNumbers.push_back(15);
	listOfSwapNumbers.push_back(25); listOfSwapNumbers.push_back(50); listOfSwapNumbers.push_back(100);
	vector<int> testSwapType;
	testSwapType.push_back(1); testSwapType.push_back(2); testSwapType.push_back(3); 
	testSwapType.push_back(4); testSwapType.push_back(5); 

	vector<int> transferInformation = getNodeLayerInfo(swapFileNames[0]);
	int numberLayers = transferInformation[0];
	int numberNodes = transferInformation[1];

	bool noMoreMotifs = false;
	int totalSwapsPerformed = 0;
	for(int k = 0; k < 30; k++){
		if(k == 6 || k==12 || k==18 || k==24) cout << endl << endl << "---------------" << endl;
		if(k == 6 || k==12 || k==18 || k==24){
			totalSwapsPerformed = 0;
			noMoreMotifs = false;
		}
		if(noMoreMotifs) continue;
		//double threshold = testingNumbers[k/6];
		swapMethod = testSwapType[k/6];
		double threshold = 0.4;
		cout << "The threshold is " << threshold << endl;
		cout << "The current swap number is " << swapMethod << endl;
		string networkName;
		if(k == 0+((k/6)*6) || k == 1+((k/6)*6)) networkName = startingNetworkName;
		else networkName = swapFileNames[k-(((k/6)*6) + 1)];
		cout << networkName << endl;

		int swapNumber = listOfSwapNumbers[k - ((k/6)*6)];

		auto f2_start = chrono::steady_clock::now();
		totalSwapsPerformed = totalSwapsPerformed + swapNumber;
		g = new Graph(numberLayers, numberNodes, networkName, threshold, swapNumber, swapMethod);
		if(g->cannotSwap){
			auto f2_end = chrono::steady_clock::now();
			//Enter failure information into the csv document
			totalCSV << swapMethod << "," << type << "," << threshold << "," << swapNumber << ",-1,-1" << endl;
			delete g;
			continue;
		}
		if(swapMethod == 5 && swapNumber > 0){
			vector<Motif*> tempPattern;
			g->findMotifs(tempPattern,type);
			string realNetworkName = g->motifShuffle(tempPattern, networkName, swapNumber);
			for(int y = 0; y < tempPattern.size(); y++) delete tempPattern[y];
			delete g;
			g = new Graph(numberLayers, numberNodes, realNetworkName, threshold, swapNumber, swapMethod);
		}

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
		totalCSV << swapMethod << "," << type << "," << threshold << "," << totalSwapsPerformed << ","<<numEmbeddings<<","<< chrono::duration_cast<chrono::milliseconds>(f2_end - f2_start).count()<< endl;

    	//save edges of final motif set
    	ofstream finalRes;
    	finalRes.open("motif_edges.txt");
    	for (int i = 0; i < selectEdge.size(); i++) {
    	    finalRes<<g->edgeset[selectEdge[i]].first<<"\t"<<g->edgeset[selectEdge[i]].second<<endl;
    	}
    	finalRes.close();
		delete g;
		for(int i = 0; i < pattern.size(); i++) delete pattern[i];
		if(numEmbeddings == 0) noMoreMotifs = true;
	}

/*
	vector<int> transferInformation = getNodeLayerInfo(networkName);
	int numberLayers = transferInformation[0];
	int numberNodes = transferInformation[1];

	auto f2_start = chrono::steady_clock::now();
	g = new Graph(numberLayers, numberNodes, networkName, threshold, swapNumber);
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
	delete g;*/

	totalCSV.close();
	
}
