from igraph import *

def lp1(max_iteration, data):
    
    for v in range(max_iteration):
        for i in range(19000):
            if data[i][1] != 0 and len(data[i][2]) > 0:
                for k in range(len(data[i][2])):
                    if data[data[i][2][k][0]][1] == 0:
                        new_label = [data[i][1], data[i][2][k][1]]
                        data[data[i][2][k][0]].append(new_label) 
                        print(data[data[i][2][k][0]])
                data[i][2] = [] 
        
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


                
