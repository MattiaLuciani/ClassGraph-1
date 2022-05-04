from igraph import *
from timeit import default_timer as timer
import memory_profiler as mem_profile
import numpy as np
import pandas as pd


def lp1(max_iteration, data):

    print('Memory (Before): ' + str(mem_profile.memory_usage()) + 'MB' )
    start = timer()

    for v in range(max_iteration): #for i = 1,2,...,#iter do
        tolabel_count = 0

        # Nodes at level i-1 send info to their neighbours at level i, populate level_1
        for i in range(len(data)): #for all j ∈ Level(i−1) do
            if data[i][1] != 0 and len(data[i][2]) > 0: #and data[i][1] > 1: #All labelled nodes with at least one egde
                for k in range(len(data[i][2])): #for all edges e_jk do
                    # if the neighbour to be labelled don't have already a label
                    if data[data[i][2][k][0]][1] == 0: #ACCORPARE #Send to the node k the label lj with weight wjk
                        tolabel_count += 1
                        new_label = [data[i][1], data[i][2][k][1]] #(l,w)
                        #new_label = [data[i][2][k][1], data[i][1]] #(w,l)
                        data[data[i][2][k][0]].append(new_label) # Add node k to Level_i
                data[i][2] = [] #Delete all the edges of node j

        # For each node at level i the final label is decided
        print(tolabel_count)

        if tolabel_count == 0:
            break

    for i in range(len(data)):
        len_line = len(data[i])
        if len_line > 3: 
            for k in range(3, len_line):
                summing_list = []
            for j in range(len(data[i][k])):
                if len(summing_list) == 0:
                    summing_list.append(data[i][k][j])
                else:
                    summing_list[len(summing_list) - 1][1] = summing_list[len(summing_list) - 1][1] + possible_labels[j][1]
            data[i][1] = summing_list[len(summing_list) - 1][0]
            for k in range(3, len_line):
                del data[i][3]

    end = timer()
    print("Time taken:", end-start)
    print('Memory (After) : ' + str(mem_profile.memory_usage()) + 'MB')


                
