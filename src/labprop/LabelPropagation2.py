from igraph import *
from timeit import default_timer as timer

def lp1(max_iteration, data):
    start = timer()
    
    for v in range(max_iteration):
        tolabel_count = 0
        for i in range(len(data)):
            if data[i][1] != 0 and len(data[i][2]) > 0:
                for k in range(len(data[i][2])):
                    if data[data[i][2][k][0]][1] == 0:
                        tolabel_count += 1
                        new_label = [data[i][1], data[i][2][k][1]]
                        data[data[i][2][k][0]].append(new_label)
                data[i][2] = []
        #print(tolabel_count)


        for i in range(len(data)):
                if len(data[i]) > 3:
                    possible_labels = []
                    for k in range(3, len(data[i])):
                        possible_labels.append(data[i][k])
                    summing_list = []
                    for j in range(len(possible_labels)):
                        if len(summing_list) == 0:
                            summing_list.append(possible_labels[j])
                        else:
                            if summing_list[len(summing_list) - 1][0] == possible_labels[j][0]:
                                summing_list[len(summing_list) - 1][1] = summing_list[len(summing_list) - 1][1] + possible_labels[j][1]
                    data[i][1] = summing_list[len(summing_list) - 1][0]
                    for k in range(3, len(data[i])):
                        del data[i][3]
    end = timer()
    print("Time taken:", end-start)
    
    
