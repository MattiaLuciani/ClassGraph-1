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

#.asqg file path
sgafile =

#Kraken2 output path
kraken2_file = 

//Kraken2 output path
output_path = 

max_iteration = 20
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

    for line in infile1:
        if line.startswith("ED"):
            words = line.split()

	    #first node
            firststring = words[1].split(".", 1)[1]

	    #second node
            secondstring = words[2].split(".", 1)[1]

	    #overlap of reads / chaining scoore
            overlap = int(words[9])
            normalizedoverlaplength = float(overlap/3000)
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
    reads_info.append(krakenlabels[i][1])
    read_id = krakenlabels[i][0]
    my_map.append(read_id)
    node_count += 1


for j in range(len(edges)):
    link = []
    first_node = edges[j][0]
    second_node = edges[j][1]
    link.append(int(first_node))
    link.append(int(second_node))
    link.append(float(edges[j][2]))
    links.append(list(link))



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
print('ok1')

# Iterate links
for i in range(len(links)):
    # Remove self loops
    if links[i][0] != links[i][1]:
        edge_list.append((int(links[i][0]), int(links[i][1])))
        
        weights_list.append((float(links[i][2])))


print("Total number of edges in the assembly graph: " + str(len(edge_list)))

label_neighbors = []

f = 0
for read in range(node_count):
    if reads_info[read] == '0': # conta le label messe a 0 
        f = f + 1 
print('label = 0 prima di RL:', f)  


#RefineLabels algorithm
to_elim = []
count = 0
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


# Run label propagation
# -----------------------

print("Preparing data for label propagation")

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

    line.append(neighs)s
    data.append(line)

print("Starting label propagation")

iteration = 0
for v in range(max_iteration):
    iteration = iteration + 1
    tolabel_count = 0
    for i in range(len(data)):
        if int(data[i][1] != 0) and len(data[i][2]) > 0:
            for k in range(len(data[i][2])):
                to_label = int(data[i][2][k][0])
                if data[to_label][1] == 0:
                    tolabel_count += 1
                    tmp = []
                    tl_weight = data[i][2][k][1]
                    label = data[i][1]
                    tmp.append(tl_weight)
                    tmp.append(label)
                    data[to_label].append(tmp)
            data[i][2] = []

    if tolabel_count == 0:
        break

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


print("***************Label propagation termined**************")

k = 0
for read in range(node_count):
    if reads_info[read] == '0': # conta le label messe a 0 
        k = k + 1 
print('label = 0 dopo di LP:', k)

output_file = output_path + 'ClassGraph2-out.res'

with open(output_file, mode='w') as out_file:
    for i in range(len(data)):
        out_file.write(krakenlabels[i][0] + "\t" + str(data[i][1]) + "\n")
        
print("Final binning results can be found at " + output_file)
