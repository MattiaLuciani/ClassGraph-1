from igraph import *
from timeit import default_timer as timer
import memory_profiler as mem_profile


def lp1(max_iteration, data):
    
    print('Memory (Before): ' + str(mem_profile.memory_usage()) + 'MB' )
    start = timer()
    
    for v in range(max_iteration): #for i = 1,2,...,#iter do

        # Nodes at level i-1 send info to their neighbours at level i, populate level_1

        for i in range(len(data)): #for all j ∈ Level(i−1) do

            if data[i][1] != 0 and len(data[i][2]) > 0: #and data[i][1] > 1: #All labelled nodes with at least one egde
                
                for k in range(len(data[i][2])): #for all edges e_jk do
                    
                    # if the neighbour to be labelled don't have already a label
                    
                    to_label = data[i][2][k][0] #Send to the node k the label lj with weight wjk
                    
                    if data[to_label][1] == 0:
                        label = data[i][0] #error data[i][1]
                        weight = data[i][2][k][1]
                        tmp = [label, weight] 
                        data[to_label].append(tmp) # Add node k to Level_i

                data[i][2] = [] #Delete all the edges of node j

        # For each node at level i the final label is decided

        for i in range(len(data)):
            
            len_line = len(data[i])
            
            if len_line > 3:
                possible_labels = []
                for k in range(3, len_line):
                    possible_labels.append(data[i][k])
                    
                summing_list = []
                
                for j in range(len(possible_labels)):
                    summing_list.append(possible_labels[j])
                    
                    #print(len(summing_list) - 1][1])

                    if summing_list[len(summing_list) - 1][1] == possible_labels[j][1]:
                        summing_list[len(summing_list) - 1][0] = summing_list[len(summing_list) - 1][0] + possible_labels[j][0]
                        
                    else:
                        summing_list.append(possible_labels[j])
                  
                data[i][1] = summing_list[len(summing_list) - 1][1] #Level_i−1 = Level_i
                
                for k in range(3, len_line): #Level_i = ∅
                    del data[i][3]
        
    end = timer()
    print("Time taken:", end-start)
    print('Memory (After) : ' + str(mem_profile.memory_usage()) + 'MB')
