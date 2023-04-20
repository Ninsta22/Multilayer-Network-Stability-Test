# Co-existingMotifsInMultilayerNetworks
This is the repository for the following paper : https://dl.acm.org/doi/abs/10.1145/3535508.3545528
1. How to run the code

 i). install c++ library 'boost'
 ii). compile the file using command 'make'
 iii). run the program using command in the format as:

    './motif Type_of_motif_pattern Network_name Threshold Swaps'
    example: ./motif 2 test.txt 1
    
    Each parameters is explained as follows:
    a.Type of motif pattern: 1 (Feedforward); 2 (Bifan); 3 (Biparallel); 4 (Cascade)
      We currently only support 4 basic motif patterns. To support large motif pattern, you can refer to paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5053092/. 

    b. Network name: path to the network file
    c. Edge Swap method: 1 (Purely Random); 2 (Non-Aggregate Graph x Non-Aggregate Graph); 3 (Aggregate Graph x Aggregate Graph); 4 (Aggregate Graph x Non-Aggregate Graph); 5 (Stochastic Shuffle)

2. Input network format
Each row is an edge: source_node tab target node
We use '-----' as a line to seperate each layer; and also add '----' at the end of file to indicate this layer of network is finished.

See test.txt file as an example.

3. output of program
   
   Given a motif type, we output 1) the number of motifs we have found, 2) the running time, 3) the layers which have these motifs, 4) a file save the edges of all these motifs (If the motif pattern has n edges, then every n edges constitue a motif embedding.
