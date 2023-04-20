#include "Graph.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include <iostream>
#include <random>
#include <bits/stdc++.h>

using namespace std;

string Graph::aggNonAggregateShuffleEdge(string networkName, int swapNumber){

	srand(time(0));
	string fileName = "";

	vector<vector<vector<int> > > nodePairTracker(LayerSize, vector<vector<int>>(NodeSize));
	unordered_map<string, dynamic_bitset<>> edgeCount;

	//Keep track of swappable non-aggregate edges
	vector<vector<pair<int, int>>> swappableNonEdges(LayerSize);

	//Keep track of swappable aggregate edges
	vector<vector<pair<int, int>>> swappableAggEdges(LayerSize);

	string str;
	ifstream in_stream;
	in_stream.open(networkName);

	int layerIndex = 0;
	while(!in_stream.eof()){
		getline(in_stream, str);
		if(str.size()>0){
			if(str[0] == '-') layerIndex++;
			else{
				vector<string> tokenizer;
				boost:split(tokenizer,str,boost::is_any_of("\t"));

				nodePairTracker[layerIndex][stoi(tokenizer[0])].push_back(stoi(tokenizer[1]));

				if(edgeCount.find(str) == edgeCount.end()){
					edgeCount.insert(make_pair(str, dynamic_bitset<>(LayerSize)));
				}
				edgeCount[str][layerIndex] = 1;
			}
		}
	}
	in_stream.close();

	vector<string> strs;
	for(auto& p: edgeCount){
		dynamic_bitset<> tdb = p.second;
		if (tdb.count() >= k) {
			string edge = p.first;
			strs.clear();
			boost::split(strs,edge,boost::is_any_of("\t"));

			int source = stoi(strs[0]);
	        int dest = stoi(strs[1]);

			for(int i = 0; i < tdb.size(); i++){
				if(tdb[i] == 1) swappableAggEdges[i].push_back(make_pair(source, dest));
			}
		}else if(tdb.count() > 0 && tdb.count() < k){
			string edge = p.first;
			strs.clear();
			boost::split(strs,edge,boost::is_any_of("\t"));

			int source = stoi(strs[0]);
			int dest = stoi(strs[1]);

			for(int i = 0; i < tdb.size(); i++){
				if(tdb[i] == 1) swappableNonEdges[i].push_back(make_pair(source, dest));
			}
		}
	}

	for(int j = 0; j < swapNumber; j++){
		bool flag = true;
		int aggIndex; int nonAggIndex;
		int randomLayer = rand() % LayerSize;
		while(swappableNonEdges[randomLayer].size() < 2) randomLayer = rand() % LayerSize;

		while(flag){
			aggIndex = rand() % swappableAggEdges[randomLayer].size();
			nonAggIndex = rand() % swappableNonEdges[randomLayer].size();

			while(swappableAggEdges[randomLayer][aggIndex].first == swappableNonEdges[randomLayer][nonAggIndex].first ||
			swappableAggEdges[randomLayer][aggIndex].second == swappableNonEdges[randomLayer][nonAggIndex].second){
				aggIndex = rand() % swappableAggEdges[randomLayer].size();
				nonAggIndex = rand() % swappableNonEdges[randomLayer].size();
			}
			flag = false;
			int shortIndex = swappableAggEdges[randomLayer][aggIndex].first;
			for(int i = 0; i < nodePairTracker[randomLayer][shortIndex].size(); i++){
				if(nodePairTracker[randomLayer][shortIndex][i] == swappableNonEdges[randomLayer][nonAggIndex].second){
					flag = true;
				}
			}
			shortIndex = swappableNonEdges[randomLayer][nonAggIndex].first;
			for(int i = 0; i < nodePairTracker[randomLayer][shortIndex].size(); i++){
				if(nodePairTracker[randomLayer][shortIndex][i] == swappableAggEdges[randomLayer][aggIndex].second){
					flag = true;
				}
			}
		}

		int node1 = swappableAggEdges[randomLayer][aggIndex].first;
		int node2 = swappableAggEdges[randomLayer][aggIndex].second;
		int node3 = swappableNonEdges[randomLayer][nonAggIndex].first;
		int node4 = swappableNonEdges[randomLayer][nonAggIndex].second;

		for(int i = 0; i < nodePairTracker[randomLayer][node1].size(); i++){
			if(nodePairTracker[randomLayer][node1][i] == node2) nodePairTracker[randomLayer][node1][i] = node4;
		}
		for(int i = 0; i < nodePairTracker[randomLayer][node3].size(); i++){
			if(nodePairTracker[randomLayer][node3][i] == node4) nodePairTracker[randomLayer][node3][i] = node2;
		}

		swappableAggEdges[randomLayer][aggIndex].second = node4;
		swappableNonEdges[randomLayer][nonAggIndex].second = node2;
	}

	fileName = to_string(swapNumber)+"swapNetwork.txt";
	ofstream swappedNetwork;
	swappedNetwork.open(fileName);

	for(auto& l: nodePairTracker){
		for(int i = 0; i < l.size(); i++){
			for(auto& p: l[i]){
				swappedNetwork << i << "\t" << p << endl;
			}
		}
		swappedNetwork << "-----------" << endl;
	}
	swappedNetwork.close();

	//Clear to prevent any leakage
	edgeCount.clear();
	nodePairTracker.clear();
	swappableAggEdges.clear();
	swappableNonEdges.clear();

	return fileName;
}

//This method will shuffle the edges corresponding to whether or not it is apart of an important
//motif
//TODO Make this return a string, which is the name of the new network
string Graph::motifShuffle(vector<Motif *> &motifs, string networkName, int swapNumber){
	srandom(time(0));
	srand(time(0));
	string fileName = "";

	vector<vector<vector<int> > > nodePairTracker(LayerSize, vector<vector<int>>(NodeSize));
	unordered_map<string, dynamic_bitset<>> edgeCount;

	//Keep track of swappable non-aggregate edges
	vector<vector<pair<int, int>>> swappableEdges(LayerSize);

	string str;
	ifstream in_stream;
	in_stream.open(networkName);

	int layerIndex = 0;
	while(!in_stream.eof()){
		getline(in_stream, str);
		if(str.size()>0){
			if(str[0] == '-') layerIndex++;
			else{
				vector<string> tokenizer;
				boost:split(tokenizer,str,boost::is_any_of("\t"));

				nodePairTracker[layerIndex][stoi(tokenizer[0])].push_back(stoi(tokenizer[1]));

				if(edgeCount.find(str) == edgeCount.end()){
					edgeCount.insert(make_pair(str, dynamic_bitset<>(LayerSize)));
				}
				edgeCount[str][layerIndex] = 1;
			}
		}
	}
	in_stream.close();

	vector<string> strs;
	for(auto& p: edgeCount){
		dynamic_bitset<> tdb = p.second;
		if(tdb.count() < k && tdb.count() > 0){
			string edge = p.first;
			strs.clear();
			boost::split(strs,edge,boost::is_any_of("\t"));

			int source = stoi(strs[0]);
			int dest = stoi(strs[1]);

			for(int i = 0; i < tdb.size(); i++){
				if(tdb[i] == 1) swappableEdges[i].push_back(make_pair(source, dest));
			}
		}
	}

	//vector<pair<int, int>> initialMotifData;
	unordered_map<int, int> initialMotifData;
	vector<pair<double, int>> motifSwapList;


	for(auto& p: motifs){
		for(int i = 0; i < p->edges.size(); i++){
			if(p->edges[i]==1){
				if(initialMotifData.find(i) == initialMotifData.end()){
					initialMotifData.insert(make_pair(i, 1));
				}else{
					initialMotifData[i] = initialMotifData[i]+1;
				}
			}
		}
	}

	double calculateImpact = 0.0;
	for(auto& p: initialMotifData){
		dynamic_bitset<> tempLayers(LayerSize);
		for(auto& t: motifs) if(t->edges[p.first] == 1) tempLayers |= t->layers;
		int numLayersPresent = tempLayers.count();
		double impactRatio = ((double)p.second/(double)numLayersPresent);
		motifSwapList.push_back(make_pair(impactRatio, p.first));
		calculateImpact = calculateImpact + impactRatio;
	}

	sort(motifSwapList.begin(), motifSwapList.end());
	int totalImpact = (int)calculateImpact;

	for(int nonUsed = 0; nonUsed < swapNumber; nonUsed++){
		int motifEdgeID = -1;
		
		bool notFound = true;
		while(notFound){
			const long max_rand = 1000000L;
			double randomPos = ((double)totalImpact)* (random() % max_rand) / max_rand;
			double target = 0.0;
			for(int x = 0; x < motifSwapList.size(); x++){
				if(target + motifSwapList[x].first < randomPos){
					target = target + motifSwapList[x].first;
				}else{
					if(edgeLayers[motifSwapList[x].second].count() == 0){
						break;
					}else{
						notFound = false;
						motifEdgeID = motifSwapList[x].second;
						break;
					}
				}
			}
		}
		
		int node1 = stoi(edgeset[motifEdgeID].first);
		int node2 = stoi(edgeset[motifEdgeID].second);
		
		
		int workingLayer = -1;
		int testingLayerValue = 0;
		bool flag = true;
		while(testingLayerValue == 0){
			workingLayer = rand() % edgeLayers[motifEdgeID].size();
			testingLayerValue = (edgeLayers[motifEdgeID][workingLayer])*(swappableEdges[workingLayer].size());
		}


		flag = true;
		int edgeIndex2 = -1;
		while(flag){
			edgeIndex2 = rand() % swappableEdges[workingLayer].size();

			//cout << nonUsed << endl;
			while(swappableEdges[workingLayer][edgeIndex2].first == node1 ||
			swappableEdges[workingLayer][edgeIndex2].second == node2){
				edgeIndex2 = rand() % swappableEdges[workingLayer].size();
			}

			flag = false;
			for(int i = 0; i < nodePairTracker[workingLayer][node1].size(); i++){
				if(nodePairTracker[workingLayer][node1][i] == swappableEdges[workingLayer][edgeIndex2].second){
					flag = true;
				}
			}
			int shortIndex = swappableEdges[workingLayer][edgeIndex2].first;
			for(int i = 0; i < nodePairTracker[workingLayer][shortIndex].size(); i++){
				if(nodePairTracker[workingLayer][shortIndex][i] == node2){
					flag = true;
					//cout << "Here second: " << node1 << " and " << node2 << endl;
				}
			}
		}

		int node3 = swappableEdges[workingLayer][edgeIndex2].first;
		int node4 = swappableEdges[workingLayer][edgeIndex2].second;

		for(int i = 0; i < nodePairTracker[workingLayer][node1].size(); i++){
			if(nodePairTracker[workingLayer][node1][i] == node2) nodePairTracker[workingLayer][node1][i] = node4;
		}
		for(int i = 0; i < nodePairTracker[workingLayer][node3].size(); i++){
			if(nodePairTracker[workingLayer][node3][i] == node4) nodePairTracker[workingLayer][node3][i] = node2;
		}

		//Make corrections to any data structures that need corrections done
		swappableEdges[workingLayer][edgeIndex2].second = node2;
		for(int i = 0; i < edgeLayers[motifEdgeID].size(); i++){
			if(edgeLayers[motifEdgeID][i] == workingLayer) edgeLayers[motifEdgeID][i] = 0;
		}
	}

	fileName = to_string(swapNumber)+"swapNetwork.txt";
	ofstream swappedNetwork;
	swappedNetwork.open(fileName);

	for(auto& l: nodePairTracker){
		for(int i = 0; i < l.size(); i++){
			for(auto& p: l[i]){
				swappedNetwork << i << "\t" << p << endl;
			}
		}
		swappedNetwork << "-----------" << endl;
	}
	swappedNetwork.close();

	//Clear to prevent any leakage
	edgeCount.clear();
	nodePairTracker.clear();
	swappableEdges.clear();
	initialMotifData.clear();
	motifSwapList.clear();

	return fileName;	
}

string Graph::aggregateShuffleEdge(string networkName, int swapNumber){

	srand(time(0));
	string fileName = "";

	vector<vector<vector<int> > > nodePairTracker(LayerSize, vector<vector<int>>(NodeSize));
	unordered_map<string, dynamic_bitset<>> edgeCount;

	//Keep track of swappable non-aggregate edges
	vector<vector<pair<int, int>>> swappableEdges(LayerSize);

	string str;
	ifstream in_stream;
	in_stream.open(networkName);

	int layerIndex = 0;
	while(!in_stream.eof()){
		getline(in_stream, str);
		if(str.size()>0){
			if(str[0] == '-') layerIndex++;
			else{
				vector<string> tokenizer;
				boost:split(tokenizer,str,boost::is_any_of("\t"));

				nodePairTracker[layerIndex][stoi(tokenizer[0])].push_back(stoi(tokenizer[1]));

				if(edgeCount.find(str) == edgeCount.end()){
					edgeCount.insert(make_pair(str, dynamic_bitset<>(LayerSize)));
				}
				edgeCount[str][layerIndex] = 1;
			}
		}
	}
	in_stream.close();

	vector<string> strs;
	for(auto& p: edgeCount) {
		dynamic_bitset<> tdb = p.second;
		if (tdb.count() >= k) {
			string edge = p.first;
			strs.clear();
			boost::split(strs,edge,boost::is_any_of("\t"));

			int source = stoi(strs[0]);
	        int dest = stoi(strs[1]);

			for(int i = 0; i < tdb.size(); i++){
				if(tdb[i] == 1) swappableEdges[i].push_back(make_pair(source, dest));
			}
		}
	}

	for(int j = 0; j < swapNumber; j++){
		bool flag = true;
		int edgeIndex1; int edgeIndex2;
		int randomLayer = rand() % LayerSize;
		while(swappableEdges[randomLayer].size() < 2) randomLayer = rand() % LayerSize;

		while(flag){
			edgeIndex1 = rand() % swappableEdges[randomLayer].size();
			edgeIndex2 = rand() % swappableEdges[randomLayer].size();

			while(edgeIndex1 == edgeIndex2 || 
			swappableEdges[randomLayer][edgeIndex1].first == swappableEdges[randomLayer][edgeIndex2].first ||
			swappableEdges[randomLayer][edgeIndex1].second == swappableEdges[randomLayer][edgeIndex2].second){
				edgeIndex2 = rand() % swappableEdges[randomLayer].size();
			}
			flag = false;
			int shortIndex = swappableEdges[randomLayer][edgeIndex1].first;
			for(int i = 0; i < nodePairTracker[randomLayer][shortIndex].size(); i++){
				if(nodePairTracker[randomLayer][shortIndex][i] == swappableEdges[randomLayer][edgeIndex2].second){
					flag = true;
				}
			}
			shortIndex = swappableEdges[randomLayer][edgeIndex2].first;
			for(int i = 0; i < nodePairTracker[randomLayer][shortIndex].size(); i++){
				if(nodePairTracker[randomLayer][shortIndex][i] == swappableEdges[randomLayer][edgeIndex1].second){
					flag = true;
				}
			}
		}

		int node1 = swappableEdges[randomLayer][edgeIndex1].first;
		int node2 = swappableEdges[randomLayer][edgeIndex1].second;
		int node3 = swappableEdges[randomLayer][edgeIndex2].first;
		int node4 = swappableEdges[randomLayer][edgeIndex2].second;

		for(int i = 0; i < nodePairTracker[randomLayer][node1].size(); i++){
			if(nodePairTracker[randomLayer][node1][i] == node2) nodePairTracker[randomLayer][node1][i] = node4;
		}
		for(int i = 0; i < nodePairTracker[randomLayer][node3].size(); i++){
			if(nodePairTracker[randomLayer][node3][i] == node4) nodePairTracker[randomLayer][node3][i] = node2;
		}

		swappableEdges[randomLayer][edgeIndex1].second = node4;
		swappableEdges[randomLayer][edgeIndex2].second = node2;
	}

	fileName = to_string(swapNumber)+"swapNetwork.txt";
	ofstream swappedNetwork;
	swappedNetwork.open(fileName);

	for(auto& l: nodePairTracker){
		for(int i = 0; i < l.size(); i++){
			for(auto& p: l[i]){
				swappedNetwork << i << "\t" << p << endl;
			}
		}
		swappedNetwork << "-----------" << endl;
	}
	swappedNetwork.close();

	//Clear to prevent any leakage
	edgeCount.clear();
	nodePairTracker.clear();
	swappableEdges.clear();

	return fileName;
}

string Graph::nonAggregateShuffleEdge(string networkName, int swapNumber){

	srand(time(0));
	string fileName = "";

	vector<vector<vector<int> > > nodePairTracker(LayerSize, vector<vector<int>>(NodeSize));
	unordered_map<string, dynamic_bitset<>> edgeCount;

	//Keep track of swappable non-aggregate edges
	vector<vector<pair<int, int>>> swappableEdges(LayerSize);

	string str;
	ifstream in_stream;
	in_stream.open(networkName);

	int layerIndex = 0;
	while(!in_stream.eof()){
		getline(in_stream, str);
		if(str.size()>0){
			if(str[0] == '-') layerIndex++;
			else{
				vector<string> tokenizer;
				boost:split(tokenizer,str,boost::is_any_of("\t"));

				nodePairTracker[layerIndex][stoi(tokenizer[0])].push_back(stoi(tokenizer[1]));

				if(edgeCount.find(str) == edgeCount.end()){
					edgeCount.insert(make_pair(str, dynamic_bitset<>(LayerSize)));
				}
				edgeCount[str][layerIndex] = 1;
			}
		}
	}
	in_stream.close();

	vector<string> strs;
	for(auto& p: edgeCount){
		dynamic_bitset<> tdb = p.second;
		if(tdb.count() < k && tdb.count() > 0){
			string edge = p.first;
			strs.clear();
			boost::split(strs,edge,boost::is_any_of("\t"));

			int source = stoi(strs[0]);
			int dest = stoi(strs[1]);

			for(int i = 0; i < tdb.size(); i++){
				if(tdb[i] == 1) swappableEdges[i].push_back(make_pair(source, dest));
			}
		}
	}

	int badLayer = 0;
	for(auto& p: swappableEdges){
		if(p.size() < 2){
			badLayer++;
		}
	}

	if(badLayer == swappableEdges.size()) cout << "CANNOT PERFORM SWAPS" << endl;
	else{
		//The following below is the actual edge swapping procedure
		//If it stalls forever, check that there are 2 unique edges to even be able to swap at the selected layer
		for(int j = 0; j < swapNumber; j++){
			bool flag = true;
			int edgeIndex1;
			int edgeIndex2;
			int randomLayer = rand() % LayerSize;
			while(swappableEdges[randomLayer].size() < 2) randomLayer = rand() % LayerSize;
			
			while(flag){
				edgeIndex1 = rand() % swappableEdges[randomLayer].size();
				edgeIndex2 = rand() % swappableEdges[randomLayer].size();

				while(edgeIndex1 == edgeIndex2 || 
				swappableEdges[randomLayer][edgeIndex1].first == swappableEdges[randomLayer][edgeIndex2].first ||
				swappableEdges[randomLayer][edgeIndex1].second == swappableEdges[randomLayer][edgeIndex2].second){
					edgeIndex2 = rand() % swappableEdges[randomLayer].size();
				}
				flag = false;
				int shortIndex = swappableEdges[randomLayer][edgeIndex1].first;
				for(int i = 0; i < nodePairTracker[randomLayer][shortIndex].size(); i++){
					if(nodePairTracker[randomLayer][shortIndex][i] == swappableEdges[randomLayer][edgeIndex2].second){
						flag = true;
					}
				}
				shortIndex = swappableEdges[randomLayer][edgeIndex2].first;
				for(int i = 0; i < nodePairTracker[randomLayer][shortIndex].size(); i++){
					if(nodePairTracker[randomLayer][shortIndex][i] == swappableEdges[randomLayer][edgeIndex1].second){
						flag = true;
					}
				}
			}

			int node1 = swappableEdges[randomLayer][edgeIndex1].first;
			int node2 = swappableEdges[randomLayer][edgeIndex1].second;
			int node3 = swappableEdges[randomLayer][edgeIndex2].first;
			int node4 = swappableEdges[randomLayer][edgeIndex2].second;

			for(int i = 0; i < nodePairTracker[randomLayer][node1].size(); i++){
				if(nodePairTracker[randomLayer][node1][i] == node2) nodePairTracker[randomLayer][node1][i] = node4;
			}
			for(int i = 0; i < nodePairTracker[randomLayer][node3].size(); i++){
				if(nodePairTracker[randomLayer][node3][i] == node4) nodePairTracker[randomLayer][node3][i] = node2;
			}

			swappableEdges[randomLayer][edgeIndex1].second = node4;
			swappableEdges[randomLayer][edgeIndex2].second = node2;
		}

		fileName = to_string(swapNumber)+"swapNetwork.txt";
		ofstream swappedNetwork;
		swappedNetwork.open(fileName);

		for(auto& l: nodePairTracker){
			for(int i = 0; i < l.size(); i++){
				for(auto& p: l[i]){
					swappedNetwork << i << "\t" << p << endl;
				}
			}
			swappedNetwork << "-----------" << endl;
		}
		swappedNetwork.close();
	}
	

	//Clear to prevent any leakage
	edgeCount.clear();
	nodePairTracker.clear();
	swappableEdges.clear();

	return fileName;
}

string Graph::randomShuffleEdge(string networkName, int swapNumber){
	
	srand(time(0));

	//vector<vector<int> > nodePairTracker(NodeSize);
	vector<vector<vector<int> > > nodePairTracker(LayerSize, vector<vector<int>>(NodeSize));
	unordered_map<string, dynamic_bitset<>> edgeCount;

	string str;
	ifstream in_stream;
    in_stream.open(networkName);


	int layerIndex = 0;
	while(!in_stream.eof()) {
		getline(in_stream, str);
		if(str.size()>0){
			if (str[0] == '-') {
				layerIndex++;
			}
			else {
				vector<string> tokenizer;
				boost::split(tokenizer,str,boost::is_any_of("\t"));

				nodePairTracker[layerIndex][stoi(tokenizer[0])].push_back(stoi(tokenizer[1]));

				if (edgeCount.find(str) == edgeCount.end()) {
					edgeCount.insert(make_pair(str, dynamic_bitset<>(LayerSize)));
				}
				edgeCount[str][layerIndex] = 1;
			}
		}
	}
	in_stream.close();
	

	//Make sure to add a user input of how many swaps to do in this for loop
	for(int j= 0; j < swapNumber; j++){
		bool flag = true;
		vector<int> sourceNodes;
		vector<int> destinationNodes;
		int randomLayer = rand()%LayerSize;
		
		while(flag){
			int firstS = rand()%NodeSize;
			int firstD = -1;
			if(nodePairTracker[randomLayer][firstS].size() > 0) firstD = nodePairTracker[randomLayer][firstS][0];

			if(firstD == -1) continue;
			else{
				sourceNodes.push_back(firstS);
				destinationNodes.push_back(firstD);
				flag = false;
			}
		}

		flag = true;
		while(flag){
			int secondS = rand()%NodeSize;
			while(secondS == sourceNodes[0] || secondS == destinationNodes[0]) secondS = rand()%NodeSize;
			bool conflict = false;
			if(nodePairTracker[randomLayer][secondS].size() > 0){
				for(int i = 0; i < nodePairTracker[randomLayer][secondS].size(); i++){
					if(nodePairTracker[randomLayer][secondS][i] == destinationNodes[0]) conflict = true;
				}
				if(!conflict){
					int randomIndex = rand()%nodePairTracker[randomLayer][secondS].size();
					sourceNodes.push_back(secondS);
					destinationNodes.push_back(nodePairTracker[randomLayer][secondS][randomIndex]);
					flag = false;
				}
			}
		}


		for(int i = 0; i < nodePairTracker[randomLayer][sourceNodes[0]].size(); i++){
			if(nodePairTracker[randomLayer][sourceNodes[0]][i] == destinationNodes[0]) nodePairTracker[randomLayer][sourceNodes[0]][i] = destinationNodes[1];
		}
		for(int i = 0; i < nodePairTracker[randomLayer][sourceNodes[1]].size(); i++){
			if(nodePairTracker[randomLayer][sourceNodes[1]][i] == destinationNodes[1]) nodePairTracker[randomLayer][sourceNodes[1]][i] = destinationNodes[0];
		}
	}

	string fileName = to_string(swapNumber)+"swapNetwork.txt";
	ofstream swappedNetwork;
	swappedNetwork.open(fileName);

	for(auto& l: nodePairTracker){
		for(int i = 0; i < l.size(); i++){
			for(auto& p: l[i]){
				swappedNetwork << i << "\t" << p << endl;
			}
		}
		swappedNetwork << "-----------" << endl;
	}

	/*for(int i = 0; i < nodePairTracker.size(); i++){
		for(auto& p: nodePairTracker[randomLayer][i]){
			swappedNetwork << i << "\t" << p << endl;
		}
	}
	swappedNetwork<<"-----------"<<endl;*/
	swappedNetwork.close();

	//Clear here to prevent leakage
    edgeCount.clear();
	nodePairTracker.clear();

	return fileName;
}

void Graph::readNetworkFile(string networkName){
	string str;
	ifstream in_stream;
	unordered_map<string, dynamic_bitset<>> edgeCount;
    in_stream.open(networkName);

//	cout <<"read file "<<fileName<<endl;
	int layerIndex = 0;
	while(!in_stream.eof()) {
		getline(in_stream, str);
		if(str.size()>0){
			if (str[0] == '-') {
				layerIndex++;
			}
			else {
				if (edgeCount.find(str) == edgeCount.end()) {
					edgeCount.insert(make_pair(str, dynamic_bitset<>(LayerSize)));
				}
				edgeCount[str][layerIndex] = 1;
			}
		}
	}
	in_stream.close();

	int sourceNode, destNode;
	vector<string> strs;
	unordered_map<string, int> nodeMap; //input string vs map id

	int edgeIndex = 0; 
	int index = 0; //for node
	for(auto& p: edgeCount) {
		dynamic_bitset<> tdb = p.second;
		if (tdb.count() >= k) {
			string edge = p.first;
			strs.clear();
			boost::split(strs,edge,boost::is_any_of("\t"));

			string s1 = strs[0];
	        string s2 = strs[1];

	        if (nodeMap.find(s1) == nodeMap.end()) {
	            nodeMap.insert(make_pair(s1, index++));
	        }

	        if (nodeMap.find(s2) == nodeMap.end()) {
	            nodeMap.insert(make_pair(s2, index++));
	        }

			sourceNode=nodeMap[s1];
			destNode=nodeMap[s2];

			if(adjacencyList[sourceNode][destNode]==0){
				adjacencyMatrix[sourceNode][destNode] = edgeIndex;
				adjacencyList[sourceNode][destNode] = 1;
				edgeLayers.push_back(tdb);
				edgeset.push_back(make_pair(s1, s2));
				edgeIntSet.push_back(make_pair(sourceNode, destNode));
//				cout<<"remaining edges "<<edgeIndex<<": "<< sourceNode<<" "<<destNode<<", "<<tdb<<endl;
				edgeIndex++;
			}
		}
	}

	EdgeSize=edgeIndex;
	strs.clear();
    edgeCount.clear();
}

void Graph::findMotifs(vector<Motif*>& pattern, int caseN){
	dynamic_bitset<>::size_type it,itt,iit;
	dynamic_bitset<> ba(LayerSize);
	dynamic_bitset<> back(LayerSize);
	dynamic_bitset<> tmp(LayerSize);

	int index1,index2,index3, index4;
	index1=index2=index3=index4=-1;

	int count = 0;
	if(caseN == 1){ //Feed forward loop
		for(int i=0;i<NodeSize;i++){
			it=adjacencyList[i].find_first();
			while(it!=dynamic_bitset<>::npos){
				index1 = adjacencyMatrix[i][it];
				ba = edgeLayers[index1];

				itt=adjacencyList[i].find_next(it);
				while(itt!=dynamic_bitset<>::npos){
					index2 = adjacencyMatrix[i][itt];
					tmp = edgeLayers[index2];
					tmp &= ba;
					if (tmp.count() >= k) {
						if(adjacencyList[it][itt]==1){
							index3=adjacencyMatrix[it][itt];
							back = edgeLayers[index3];
							back &= tmp;

							if (back.count() >= k) {
								dynamic_bitset<> edges(EdgeSize);
								edges[index1] = edges[index2] = edges[index3] = 1;
								Motif* t=new Motif(count, edges, back);
								pattern.push_back(t);
								count++;
							}
						}
						if(adjacencyList[itt][it]==1){
							index3=adjacencyMatrix[itt][it];
							back = edgeLayers[index3];
							back &= tmp;

							if (back.count() >= k) {
								dynamic_bitset<> edges(EdgeSize);
								edges[index1] = edges[index2] = edges[index3] = 1;
								Motif* t=new Motif(count, edges, back);
								pattern.push_back(t);
								count++;
							}
						}
					}
					itt=adjacencyList[i].find_next(itt);
				}
				it=adjacencyList[i].find_next(it);
			}
		}
	}
	else if(caseN==2){ //bifan
		for(int i=0;i<NodeSize-1;i++){
			if (adjacencyList[i].count() < 2) {
				continue;
			}
			for(int j=i+1;j<NodeSize;j++){
				if (adjacencyList[j].count() < 2) {
					continue;
				}
				dynamic_bitset<> tdb = adjacencyList[i];
				tdb &= adjacencyList[j];
				tdb[i] = 0;
				tdb[j] = 0;
				if(tdb.count()>=2){
//					cout << "consider nodes "<< i <<" and "<<j<<", their common neighbors: "<<tdb<<endl;
					it = tdb.find_first();
					while(it!=dynamic_bitset<>::npos){
//						cout<<"their first common node: "<<it<<endl;
						index1 = adjacencyMatrix[i][it];
						index2 = adjacencyMatrix[j][it];
						tmp = edgeLayers[index1];
						tmp &= edgeLayers[index2];
//						cout<<"edges: "<<index1 <<" and "<<index2<<": their common layers: "<<tmp<<endl;

						if (tmp.count() >= k) {
							itt=tdb.find_next(it);
							while(itt!=dynamic_bitset<>::npos){
//								cout<<"their next common node: "<<itt<<endl;
								index3 = adjacencyMatrix[i][itt];
								index4 = adjacencyMatrix[j][itt];
								back = edgeLayers[index3];
								back &= edgeLayers[index4];
//                cout<<"edges: "<<index3 <<" and "<<index4<<": their common layers: "<<back<<endl;

								if (back.count() >= k) {
									back &= tmp;
									if (back.count() >= k) {
//										cout<<"generate a motif"<<endl;
                                        dynamic_bitset<> edges(EdgeSize);
                                        edges[index1] = edges[index2] = edges[index3] = edges[index4] = 1;
										Motif* t=new Motif(count, edges, back);
										pattern.push_back(t);
//										cout<<"Motif "<<count<<":\t"<<edgeset[index1].first<<"\t"<<edgeset[index2].first<<"\t"<<edgeset[index1].second<<"\t"<<edgeset[index3].second<<", "<<back<<endl;
//										cout<<"Motif "<<count<<": "<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<"\t"<<", layers: "<< back<<endl;
										count++;
									}
								}
								else {
									tdb[itt] = 0;
								}
								itt = tdb.find_next(itt);
							}
						}
						else {
							tdb[it] = 0;
						}
						it = tdb.find_next(it);
					}
				}
			}
		}
	}
	else if(caseN==3){ //bi-parallel
		for(int i=0;i<NodeSize;i++){
			if (adjacencyList[i].count() < 2) {
				continue;
			}
//			cout<<"consider node "<< i << endl;
			it = adjacencyList[i].find_first();
			while(it!=dynamic_bitset<>::npos){
				index1=adjacencyMatrix[i][it];
				tmp = edgeLayers[index1];
				itt = adjacencyList[i].find_next(it);
				while(itt != dynamic_bitset<>::npos){
					index2=adjacencyMatrix[i][itt];
					back = edgeLayers[index2];
					back &= tmp;

					if (back.count() >= k) {
//            cout<< "consider other two nodes: "<<it<<" "<<itt<<", their common layers: "<<back<<endl;
						dynamic_bitset<> tdb=adjacencyList[it];
						tdb &= adjacencyList[itt];
						tdb[i]=0;
						tdb[it] = 0;
						tdb[itt] = 0;

						iit = tdb.find_first();
						while(iit!=dynamic_bitset<>::npos){
							index3=adjacencyMatrix[it][iit];
							index4=adjacencyMatrix[itt][iit];

							ba = edgeLayers[index3];
							ba &= edgeLayers[index4];
							ba &= back;
//							cout<<"consider fourth node: "<< iit <<", their common layers: "<< ba <<endl;
							if (ba.count() >= k) {
//								cout <<"generate a motif"<<endl;
                                dynamic_bitset<> edges(EdgeSize);
                                edges[index1] = edges[index2] = edges[index3] = edges[index4] = 1;
								Motif* t=new Motif(count, edges, ba);
								pattern.push_back(t);
//								cout<<"Motif "<<count<<":\t"<<edgeset[index1].second<<"\t"<<edgeset[index2].second<<"\t"<<edgeset[index1].first<<"\t"<<edgeset[index3].second<<", "<<ba<<endl;
//								cout<<"Motif "<<count<<": "<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<"\t"<<", layers: "<< ba<<endl;
								count++;
							}
							iit = tdb.find_next(iit);
						}
					}

					itt = adjacencyList[i].find_next(itt);
				}
				it = adjacencyList[i].find_next(it);
			}
		}
	}

	else if(caseN==4){ // cascade and delay
		for(int i=0;i < NodeSize - 2;i++){ //i is the smallest number among three nodes
			it = adjacencyList[i].find_next(i);
			while(it!=dynamic_bitset<>::npos){
				index1=adjacencyMatrix[i][it];
				tmp = edgeLayers[index1];
				itt = adjacencyList[it].find_next(i);
				while(itt != dynamic_bitset<>::npos){
					if(itt != it && adjacencyList[itt][i]==1){
						index2=adjacencyMatrix[it][itt];
						index3=adjacencyMatrix[itt][i];

						back = edgeLayers[index2];
						back &= edgeLayers[index3];
						back &= tmp;
//						cout <<"consider nodes: " <<i<<" "<<it<<" "<<itt<<", and their common layers "<<back<<endl;
						if (back.count() >= k) {
//							cout<<"generate a motif"<<endl;
                            dynamic_bitset<> edges(EdgeSize);
                            edges[index1] = edges[index2] = edges[index3] = 1;
							Motif* t=new Motif(count, edges, back);
							pattern.push_back(t);
							count++;
						}
					}
					itt = adjacencyList[it].find_next(itt);
				}
				it = adjacencyList[i].find_next(it);
			}
		}

	}
}
