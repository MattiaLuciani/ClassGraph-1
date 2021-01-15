from igraph import *


def lp1(max_iteration, data):
    for v in range(max_iteration):
        # print("iter" + str(v))

        for i in range(len(data)):
            if int(data[i][1] != 0) and len(data[i][2]) > 0:
                for k in range(len(data[i][2])):
                    # if the neighbour to be labelled don't have already a label
                    to_label = int(data[i][2][k][0])
                    if data[to_label][1] == 0:
                        tmp = []
                        tl_weight = data[i][2][k][1]
                        label = data[i][1]
                        tmp.append(tl_weight)
                        tmp.append(label)
                        data[to_label].append(tmp)
                data[i][2] = []

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


def lp2(max_iteration, data):
    for v in range(max_iteration):
        # print("iter" + str(v))
        tolabel_count = 0
        for i in range(len(data)):
            if int(data[i][1] != 0) and len(data[i][2]) > 0:
                for k in range(len(data[i][2])):
                    # if the neighbour to be labelled don't have already a label
                    to_label = int(data[i][2][k][0])
                    if data[to_label][1] == 0:
                        tolabel_count += 1
                        tmp = []
                        tl_weight = data[i][2][k][1]
                        label = data[i][1]
                        tmp.append(tl_weight)
                        tmp.append(label)
                        data[to_label].append(tmp)
                        data[to_label][3] += 1
                        # print(data[to_label][3])
                data[i][2] = []

        if tolabel_count == 0:
            break

        for i in range(len(data)):
            if int(data[i][1] == 0) and len(data[i]) > 4:
                for k in range(len(data[i][2])):
                    if data[data[i][2][k][0]][1] == 0 and len(data[data[i][2][k][0]]) > 4:
                        for j in range(4, 4 + data[i][3]):
                            new_tmp = []
                            mweight = data[i][j][0] * data[i][2][k][1] * 0.9
                            new_tmp.append(mweight)
                            new_tmp.append(data[i][j][1])
                            data[data[i][2][k][0]].append(new_tmp)

        for i in range(len(data)):
            len_line = len(data[i])
            if len_line > 4:
                possible_labels = []
                for k in range(4, len_line):
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
                data[i][3] = 0
                for k in range(4, len_line):
                    del data[i][4]
