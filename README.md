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
ClassGraph requires two input ﬁles: one representing a graph of reads and the other containing the result of the classiﬁcation process. The labelles assigned to the reads by the pre-existing binning tool will be propagated over the graph to the still unclassified reads.

### Overlap Graph
The graph must be composed as follows: the nodes are represented by the reads and the edges by their overlap. The file must be presented in asqg format.
One possible tools that satisfies this requirements is [**SGA**](https://github.com/jts/sga) (String Graph Assembler), an assembler based on the overlap-layout-consensus.
