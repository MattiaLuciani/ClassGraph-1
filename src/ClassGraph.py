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
from labprop.LabelPropagation import lp1, lp2




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

# Setup argument parser
# -----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--graph", required=True, help="path to the reads graph file")
ap.add_argument("--binned", required=True,
                help="path to the file with the initial binning output")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=True, help="prefix for the output file")
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

# Setup output path for log file
# ---------------------------------------------------

fileHandler = logging.FileHandler(output_path + "/" + prefix + "classgraph.log")
fileHandler.setLevel(logging.INFO)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)

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

with open(sgafile, 'r') as infile1, open(kraken2_file, 'r') as infile2:
    for line2 in infile2:
        krakenlabel = []
        krakenlabel.append(line2.split()[0])
        krakenlabel.append(line2.split()[1])
        # if line2.split()[1] == "1":
        #     print(line2.split()[0])
        krakenlabels.append(krakenlabel)
    for line in infile1:
        if line.startswith("VT"):
            words = line.split()
            readslen = len(words[2])
            break
    # print(readslen)
    prevsuffix = 0
    for line in infile1:
        if line.startswith("ED"):
            words = line.split()
            firststring = words[1]
            secondstring = words[2]
            firststringlen = len(firststring)
            secondstringlen = len(secondstring)
            newfirststr = firststring[:firststringlen-2]
            newsecondstr = secondstring[:secondstringlen - 2]
            overlapstart = int(words[3])
            overlapend = int(words[4])
            normalizedoverlaplength = float((overlapend - overlapstart + 1)/readslen)
            tmpsuffix = newfirststr.split(".", 1)[1]
            if tmpsuffix != prevsuffix:
                prevsuffix = tmpsuffix
                edge = []
                edge.append(newfirststr)
                edge.append(newsecondstr)
                edge.append(normalizedoverlaplength)
                edges.append(edge)
                subgroupcnt = 1
            else:
                alreadypresent = False
                for k in range(len(edges)-1-subgroupcnt, len(edges)):
                    if edges[k][1] == newsecondstr:
                        edges[k][2] = max(edges[k][2], normalizedoverlaplength)
                        alreadypresent = True
                        break
                if alreadypresent == False:
                    edge = []
                    edge.append(newfirststr)
                    edge.append(newsecondstr)
                    edge.append(normalizedoverlaplength)
                    edges.append(edge)
                    subgroupcnt += 1


# Get the links from the .asqg file
# -----------------------------------

links = []

my_map = []

node_count = 0

binned_reads = []

#binned_reads_info = []
reads_info = []

read_groups = []

try:
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

    # for i in range(len(read_groups)):
    #     print(read_groups[i])


    for j in range(len(edges)):
        link = []
        common1 = edges[j][0].split(".", 1)[0]
        for i in range(len(read_groups)):
            if common1 == read_groups[i][0]:
                link.append(int(read_groups[i][1] + int(edges[j][0][(len(common1) + 1):]) - 1))
                break
        common2 = edges[j][1].split(".", 1)[0]
        for i in range(len(read_groups)):
            if common2 == read_groups[i][0]:
                link.append(int(read_groups[i][1] + int(edges[j][1][(len(common2) + 1):]) - 1))
                break
        link.append(float(edges[j][2]))
        links.append(list(link))

# for i in range(len(links)):
#     print(links[i])

except:
    logger.error("Please make sure that the correct path to the assembly graph file is provided.")
    logger.info("Exiting ReadGraph!")
    sys.exit(1)

# Construct the assembly graph
# -------------------------------
#
try:

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

except:
    logger.error("Please make sure that the correct path to the assembly graph file is provided.")
    logger.info("Exiting ClassGraph... Bye...!")
    sys.exit(1)

logger.info("Total number of edges in the assembly graph: " + str(len(edge_list)))


# Run label propagation
# -----------------------

logger.info("Preparing data for label propagation")

# In the graph are not inserted edges that connect already labbeled nodes

data = []
LabbeledVertices = []

for read in range(node_count):

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

output_file = output_path + prefix + '.res'

# with open(output_file, mode='w') as out_file:
#     for i in range(len(data)):
#         read_seqid = data[i][0]
#         rg_len = len(read_groups)
#         if rg_len == 1:
#             read_id = read_groups[0][0] + "." + str(read_seqid + 1)
#         elif rg_len > 1:
#             for k in range(rg_len):
#                 if read_seqid >= read_groups[rg_len-k-1][1]:
#                     read_id = read_groups[rg_len-k-1][0] + "." + str(read_seqid - read_groups[rg_len-k-1][1] + 1)
#                     break
#         out_file.write(read_id + "\t" + str(data[i][1]) + "\n")

with open(output_file, mode='w') as out_file:
    for i in range(len(data)):
        out_file.write(krakenlabels[i][0] + "\t" + str(data[i][1]) + "\n")

logger.info("Final binning results can be found at " + output_file)
