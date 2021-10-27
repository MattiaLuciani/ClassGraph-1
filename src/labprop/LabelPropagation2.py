from igraph import *
from timeit import default_timer as timer

# First version of the propagation algorithm

def lp1(max_iteration, data):
    
    #file_uno = open("data.txt", "w")
    #str1 = ' '.join(map(str, data))
    #file_uno.write(str1)
    #file_uno.close()
    #print(data)
    start = timer()
    for v in range(max_iteration):
        # Nodes at level i-1 send info to their neighbours at level i
        tolabel_count = 0
        
        for i in range(len(data)):
            if int(data[i][1] != 0) and len(data[i][2]) > 0:
                for k in range(len(data[i][2])):
                    # if the neighbour to be labelled don't have already a label
                    to_label = int(data[i][2][k][0])
                    if data[to_label][1] == 0: # node to be labelled
                        label = data[i][1]
                        tl_weight = data[i][2][k][1]
                        tmp = [label, tl_weight]
                        data[to_label].append(tmp)
                        #tmp.append(label)
                        #tmp.append(tl_weight)
                        #tmp.sort(reverse=True)

                        #for i in range(len(tmp)):
                        
                        #tmp.reverse()
                        
                        #print(tmp)
                        
                        
                        #print(data[to_label])
        #print(tolabel_count)         
                #data[i][2] = []

            
            
        # For each node at level i the final label is decided

        for i in range(len(data)):
            len_line = len(data[i])
            if len_line > 3:
                possible_labels = []
                for k in range(3, len_line):
                    possible_labels.append(data[i][k])
                summing_list = []
                for j in range(len(possible_labels)):
                    #print(possible_labels[j])
                    summing_list.append(possible_labels[j])
                    if len(summing_list) > 0 and summing_list[len(summing_list) - 1][1] == possible_labels[j][1]:
                        summing_list[len(summing_list) - 1][0] = summing_list[len(summing_list) - 1][0] + possible_labels[j][0]
                    else:
                        summing_list.append(possible_labels[j])
                        #print("****************")
                data[i][1] = summing_list[len(summing_list) - 1][1]
               # print(data[i][1])
                for k in range(3, len_line):
                    del data[i][3]
    end = timer()
    print("Time taken:", end-start)
