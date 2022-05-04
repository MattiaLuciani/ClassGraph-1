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
import cairocffi as cairo
import igraph as ig
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

sgafile = '/Users/mattialuciani/Desktop/università_magistrale/bioinformatics/tesi/Sim_8/Sim8.asqg'


kraken2_file = '/Users/mattialuciani/Desktop/università_magistrale/bioinformatics/tesi/kraken2_Sim_8_DBKraken2.res'

# kraken
output_path = '/Users/mattialuciani/ClassGraph/output'
#prefix = 'centrifuge_RL2'
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

    for line in infile1:
        if line.startswith("ED"):
            words = line.split()

            firststring = words[1]
            secondstring = words[2]
            overlap = int(words[9])
            normalizedoverlaplength = float(overlap/100)
            
            tmpsuffix = firststring.split("_", 1)[1]
            if tmpsuffix != prevsuffix:
                prevsuffix = tmpsuffix
                edge = []
                edge.append(firststring)
                edge.append(secondstring)
                edge.append(normalizedoverlaplength)
                edges.append(edge)
                subgroupcnt = 1
            else:
                alreadypresent = False
                for k in range(len(edges)-1-subgroupcnt, len(edges)):
                    if edges[k][1] == secondstring:
                        edges[k][2] = max(edges[k][2], normalizedoverlaplength)
                        alreadypresent = True
                        break
                if alreadypresent == False:
                    edge = []
                    edge.append(firststring)
                    edge.append(secondstring)
                    edge.append(normalizedoverlaplength)
                    edges.append(edge)
                    subgroupcnt += 1

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

for j in range(len(edges)):
    link = []
    #vedere meglio i conrolli nella versione originale 
    first_node = edges[j][0].split("_", 1)[1]
    second_node = edges[j][1].split("_", 1)[1]
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
print(len(links))



#for i in range(len(links)):
for i in range(1965562):
    
    # Remove self loops
    if links[i][0] != links[i][1]:
        edge_list.append((int(links[i][0]), int(links[i][1])))
        weights_list.append((float(links[i][2])))
#print(edge_list)

print(links[1965562])

reads_graph.add_edges(edge_list)

reads_graph.es["weight"] = weights_list


logger.info("Total number of edges in the assembly graph: " + str(len(edge_list)))

label_neighbors = []
count = 0 
h = 0

'''
#first implementation RefineLabels

for read in range(node_count): 
    neighbours = reads_graph.neighbors(read)

    for v in neighbours:
        if reads_info[v] != reads_info[read] and reads_info[read] != '0' and reads_info[v] != '0':
            reads_info[read] = '0'
            reads_info[v] = '0'
'''

                

              
'''
#second implementation RefineLabels

for read in range(node_count):
    neighbours = reads_graph.neighbors(read)
    n = []
    
    if len(neighbours) > 2 and len(neighbours) % 2 != 0: #seleziono nodi con numero di vicini dispari
        
        for v in neighbours:
            n.append(reads_info[v])
            
        count = n.count(reads_info[read])
        if count < len(n)/2: 
            reads_info[read] = '0'
        #print(read, reads_info[read], n)

for read in range(node_count):
    if reads_info[read] == '0': # conta le label messe a 0 
        h = h + 1 
print('zeros', h)
'''










#refine_label(kraken2_file)

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
try:

    if labprop_v == 1:
        lp1(max_iteration, data)
    else:
        lp2(max_iteration, data)

except:
    logger.error("Please make sure that you inserted the correct parameter for the lp version (either 1 or 2)")
    logger.info("Exiting ClassGraph.. Bye...!")
    sys.exit(1)



logger.info("***************Label propagation termined**************")

#output_file = output_path + prefix + '.res'
output_file = output_path + 'Sim8_DBKraken2_CG.res'

with open(output_file, mode='w') as out_file:
    for i in range(len(data)):
        out_file.write(krakenlabels[i][0] + "\t" + str(data[i][1]) + "\n")
        
logger.info("Final binning results can be found at " + output_file)