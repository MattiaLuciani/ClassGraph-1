import sys
import os
import subprocess
import csv
import operator
import time
import argparse
import re
import logging
from igraph import * 
import igraph as ig
import collections
#from labprop.LabelPropagation import lp1
from labprop.LabelPropagation import lp1
#from labprop.LabelPropagation3 import lp1
#from reflab.RefineLabel import refine_label





# Sample command
# -------------------------------------------------------------------
# python ReadGraph_SGA.py     --graph /path/to/graph_file.asqg
#                            --binned /path/to/binning_result.csv
#                            --output /path/to/output_folder
# -------------------------------------------------------------------


# Setup logger
# -----------------------

logger = logging.getLogger('ClassGraph')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
logger.addHandler(consoleHeader)

start_time = time.time()
"""
# Setup argument parser
# -----------------------
ap = argparse.ArgumentParser()
ap.add_argument("--graph", required=True, help="path to the reads graph file")
ap.add_argument("--binned", required=True,
                help="path to the file with the initial binning output")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=False, help="prefix for the output file")
ap.add_argument("--max_iteration", required=False, type=int, default=20,
                help="maximum number of iterations for label propagation algorithm. [default: 20]")
ap.add_argument("--lp_version", required=False, type=int, default=1,
                help="Type 1 if you want to propagate with lp-v1, type 2 if you prefer to use lp-v2. [default 1]")
args = vars(ap.parse_args())
sgafile = args["graph"]
# asqg
kraken2_file = args["binned"]
# kraken
output_path = args["output"]
prefix = args["prefix"]
max_iteration = args["max_iteration"]
labprop_v= args["lp_version"]
"""

#sgafile = '/Users/mattialuciani/Desktop/università_magistrale/bioinformatics/tesi/Sim8/Sim8-1.asqg'


kraken2_file = '/Users/mattialuciani/Desktop/università_magistrale/bioinformatics/tesi/Sim8/kraken2_Sim_8_DBKraken2.res'
#kraken2_file = '/nfsd/bcb/bcbg/luciani/Kraken2/Results/kraken2_Sim_50_strex.res'

# kraken
output_path = '/Users/mattialuciani/Desktop/università_magistrale/bioinformatics/tesi/Sim8'
#output_path = '/nfsd/bcb/bcbg/luciani/CG_Sim_Results/'

max_iteration = 20
labprop_v= 1

# Setup output path for log file
# ---------------------------------------------------
#fileHandler = logging.FileHandler(output_path + "/" + prefix + "classgraph.log")
'''
fileHandler = logging.FileHandler(output_path + "/" + prefix)
fileHandler.setLevel(logging.INFO)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)


# Setup output path for log file
# ---------------------------------------------------
fileHandler = logging.FileHandler(output_path + "/" + prefix + "classgraph.log")
fileHandler.setLevel(logging.INFO)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)
'''
logger.info(
    "ReadGraph makes use of the assembly graph produced by SGA")

logger.info("Overlap graph file: " + sgafile)
logger.info("Existing binning output file: " + kraken2_file)
logger.info("Final binning output file: " + output_path)
logger.info("Maximum number of iterations: " + str(max_iteration))
logger.info("Label propagation v" + str(labprop_v))
logger.info("ReadGraph started")

# Get the number of bins from the initial binning result
# --------------------------------------------------------

logger.info("Combining the graph file with the classification one, in order to obtain a labelled graph")

krakenlabels = []
edges = []

lenght = [] 
with open(sgafile, 'r') as infile1, open(kraken2_file, 'r') as infile2:
    for line2 in infile2:
        krakenlabel = []
        krakenlabel.append(line2.split()[0])
        krakenlabel.append(line2.split()[1])
        krakenlabels.append(krakenlabel)

    prevsuffix = 0
    tmpsuffix = 0

#marine dataset
    for line in infile1:
        if line.startswith("ED"):
            words = line.split()
            firststring = words[1]
            secondstring = words[2]
            overlap = int(words[9])
            normalizedoverlaplength = float(overlap/300)
            tmpsuffix = firststring.split("S0R.")
            edge = []
            edge.append(firststring)
            edge.append(secondstring)
            edge.append(normalizedoverlaplength)
            edges.append(edge)


# Get the links from the .asqg file
# -----------------------------------

links = []

my_map = []

node_count = 0

binned_reads = []

reads_info = []

read_groups = []


common = ""
commonend = "."

for i in range(len(krakenlabels)):
    actual_common = krakenlabels[i][0].split(".", 1)[0]
    if common != actual_common:
        read_group = []
        read_group.append(actual_common)
        read_group.append(node_count)
        read_groups.append(read_group)
        common = actual_common

    if krakenlabels[i][1] != "0":
        binned_reads.append(node_count)

    reads_info.append(krakenlabels[i][1])
    read_id = krakenlabels[i][0]
    my_map.append(read_id)
    node_count += 1

           # tmpsuffix = firststring.split("S0R")
          #  if tmpsuffix[1] != prevsuffix:
for j in range(len(edges)):
    link = []
    first_node = edges[j][0].split("S0R")[1]
    second_node = edges[j][1].split("S0R")[1]
    link.append(int(first_node))
    link.append(int(second_node))
    link.append(float(edges[j][2]))
    links.append(list(link))
# Construct the assembly graph
# -------------------------------
#

# Create the graph
reads_graph = Graph()

# Create list of edges
edge_list = []

weights_list = []

# Add vertices
reads_graph.add_vertices(node_count)
# Name vertices
for i in range(len(reads_graph.vs)):
    reads_graph.vs[i]["id"] = i
    reads_graph.vs[i]["label"] = str(my_map[i])

# Iterate links
for i in range(len(links)):
    # Remove self loops
    if links[i][0] != links[i][1]:
        edge_list.append((int(links[i][0]), int(links[i][1])))
        weights_list.append((float(links[i][2])))
reads_graph.add_edges(edge_list)

reads_graph.es["weight"] = weights_list


logger.info("Total number of edges in the assembly graph: " + str(len(edge_list)))

label_neighbors = []

f = 0
for read in range(node_count):
    if reads_info[read] == '0': # conta le label messe a 0 
        f = f + 1 
print('label = 0 prima di RL:', f)  


#first implementation RefineLabels

for read in range(node_count): 
    neighbours = reads_graph.neighbors(read)
    for v in neighbours:
        if reads_info[v] != reads_info[read] and reads_info[read] != '0' and reads_info[v] != '0':
            reads_info[read] = '0'
            reads_info[v] = '0'
       
            



'''
#second implementation RefineLabels

for read in range(node_count):
    neighbours = reads_graph.neighbors(read)
    n = []        
    for v in neighbours:
        if reads_info[v] != '0':
            n.append(reads_info[v])
            count = n.count(reads_info[read])
            if count < len(n)/2: 
                reads_info[read] = '0' 
'''



'''



#third implementation RefineLabels
to_elim = []

for read in range(node_count):
    neighbours = reads_graph.neighbors(read)
    n = []        
    for v in neighbours: 
        if reads_info[v] != '0': 
            n.append(reads_info[v]) 
            count = n.count(reads_info[read]) 
    if count < len(n)/2:
        to_elim.append(read)

for read in range(node_count):
    if read in to_elim: 
        reads_info[read] = '0'


'''


'''
to_elim = []

for read in range(node_count):
    neighbours = reads_graph.neighbors(read)
    n = []        
    for v in neighbours: 
        if reads_info[v] != '0': 
            n.append(reads_info[v]) 
            count = n.count(reads_info[read]) 
    if count < len(n)/2:
        to_elim.append(read)
        
for read in to_elim: 
    reads_info[read] = '0'
'''
    
    
g = 0 
for read in range(node_count):
    if reads_info[read] == '0': # conta le label messe a 0 
        g = g + 1 
print('label = 0 dopo di RL:', g)

'''
for read in to_elim:  #3, 6, 8, 9, 12, 13, 15, 17, 20, 24, 25, 28, 32, 38
    reads_info[read] = '0'
    
#third implementation RefineLabels

for read in range(node_count):
    line = []
    line.append(read)
    line.append(int(reads_info[read]))
    neighbours = reads_graph.neighbors(read)
    neighs = []
    for neighbour in neighbours:
        n = []
        n.append(neighbour)
        n.append(reads_info[neighbour])
        e_weight = reads_graph.es[reads_graph.get_eid(read, neighbour)]["weight"]
        n.append(e_weight)
        neighs.append(n)
    line.append(neighs)
    
    

print(line) #[128692, 1515, [[114324, '1515', 0.96], [114860, '1515', 0.44333333333333336], [115882, '203119', 0.37333333333333335], [115955, '1515', 0.49333333333333335], [116032, '1515', 0.5033333333333333], [116709, '1515', 0.47], [116729, '9606', 0.5766666666666667], [117335, '1515', 0.5], [117456, '1515', 0.5266666666666666], [117508, '0', 0.67], [117988, '1515', 0.78], [118059, '1515', 0.56], [119446, '1515', 0.4166666666666667], [119530, '2', 0.4033333333333333], [121087, '1515', 0.4266666666666667], [122014, '1515', 0.43], [123148, '203119', 0.5066666666666667], [123234, '1515', 0.39], [123919, '1515', 0.39], [124592, '1515', 0.4766666666666667], [125143, '1515', 0.44], [125227, '1515', 0.47], [126500, '1515', 0.49666666666666665], [127078, '1515', 0.48], [127097, '203119', 0.37666666666666665]]]

'''



# Run label propagation
# -----------------------

logger.info("Preparing data for label propagation")

# In the graph are not inserted edges that connect already labbeled nodes
LabbeledVertices = []
data = []
count = 0
for read in range(node_count):
    #count = count + 1
    line = []
    line.append(read)
    
    if int(reads_info[read]) != 0: 
        alreadyLabelled = []
        alreadyLabelled.append(read)
        alreadyLabelled.append(int(reads_info[read]))
        LabbeledVertices.append(alreadyLabelled)


    line.append(int(reads_info[read]))
    neighbours = reads_graph.neighbors(read)
    neighs = []

    for neighbour in neighbours:
        if int(reads_info[neighbour]) == 0:
            
            n = []
            n.append(neighbour)
            e_weight = reads_graph.es[reads_graph.get_eid(read, neighbour)]["weight"]
            n.append(e_weight)
            neighs.append(n)

    line.append(neighs)
    if labprop_v == 2:
        line.append(0)
    data.append(line)

logger.info("Starting label propagation")

# The propagation at each iteration is performed from the last labelled nodes to their neighbors
# Once a node is labbeled , that label is not going to be changed


h = 0
for read in range(node_count):
    if reads_info[read] == '0': # conta le label messe a 0 
        h = h + 1 
print('label = 0 prima di LP:', h)

for v in range(max_iteration):

    # Nodes at level i-1 send info to their neighbours at level i
    tolabel_count = 0
    for i in range(len(data)):
        if int(data[i][1] != 0) and len(data[i][2]) > 0:
            for k in range(len(data[i][2])):
                # if the neighbour to be labelled don't have already a label
                to_label = int(data[i][2][k][0])
                if data[to_label][1] == 0:
                    print(data[to_label])

                    tolabel_count += 1
                    tmp = []
                    tl_weight = data[i][2][k][1]
                    label = data[i][1]
                    tmp.append(tl_weight)
                    tmp.append(label)
                    data[to_label].append(tmp)
            data[i][2] = []

    #print(tolabel_count)

    if tolabel_count == 0:
        break
    # For each node at level i the final label is decided

    for i in range(len(data)):
        len_line = len(data[i])
        if len_line > 3:
            possible_labels = []
            for k in range(3, len_line):
                possible_labels.append(data[i][k])

            possible_labels = sorted(possible_labels, key=operator.itemgetter(1))

            summing_list = []

            for j in range(len(possible_labels)):
                if len(summing_list) == 0:
                    summing_list.append(possible_labels[j])
                elif len(summing_list) > 0 and summing_list[len(summing_list) - 1][1] == possible_labels[j][1]:
                    summing_list[len(summing_list) - 1][0] = summing_list[len(summing_list) - 1][0] + \
                                                             possible_labels[j][0]
                else:
                    summing_list.append(possible_labels[j])
            summing_list = sorted(summing_list, key=operator.itemgetter(0))
            data[i][1] = summing_list[len(summing_list) - 1][1]
            for k in range(3, len_line):
                del data[i][3]





logger.info("***************Label propagation termined**************")

k = 0
for read in range(node_count):
    if reads_info[read] == '0': # conta le label messe a 0 
        k = k + 1 
print('label = 0 dopo di LP:', k)

#output_file = output_path + prefix + '.res'
output_file = output_path + 'finale.res'

with open(output_file, mode='w') as out_file:
    for i in range(len(data)):
        out_file.write(krakenlabels[i][0] + "\t" + str(data[i][1]) + "\n")
        
logger.info("Final binning results can be found at " + output_file)
