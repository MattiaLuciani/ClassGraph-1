# ClassGraph
ClassGraph is a tool that classifies metagenomic reads starting from the output of other pre-existing binning tools
## ClassGraph Download
It's possible to download ClassGraph by cloning the repository in your machine.

```
git clone https://github.com/MargheritaCavattoni/ClassGraph.git
```
## Dependencies
The istallation of ClassGraph requires python 3.6 or above. Besides the following dependency is needed:
* [python-igraph](https://igraph.org/python/)

## Preprocessing
ClassGraph requires two input ﬁles: one representing a graph of reads and the other containing the result of the classiﬁcation process. The labelles assigned to the reads by the pre-existing binning tool will be propagated over the graph to the still unclassified reads. ClassGraph is thought to be used with paired-end reads.

### Overlap Graph
The graph must be composed as follows: the nodes represent the reads and the edges their overlaps. The file must be presented in asqg format.
One possible tool that satisfies these requirements is [**SGA**](https://github.com/jts/sga) (String Graph Assembler), an assembler based on the overlap-layout-consensus.

### Binning Tools
ClassGraph requires as input the result of the classification process presented as follows:
```
NODE_1,L1
NODE_2,L1
NODE_3,L1
NODE_4,L2
NODE_5,L2
...
```
Where NODE_i is the node ID and Li represents its label. Li must be equal to 0 if the read wasn't classified.
Hypothetically any metagenomic binning tool could be used for the classification. During the testing we decided to use [**kraken2**](https://github.com/DerrickWood/kraken2.git), since it's one of the best performing.

## Usage Example
```
python3 ReadGraph.py --graph $SGA_DIR/Graph.asqg --output $Read-Graph_Output_DIR/ --binned $Kraken2_DIR/BinnedReads.out --preﬁx example1 --max_iteration 20
```
