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


def refine_label(kraken2_file):

    closest_neighbours = assembly_graph.neighbors(i, mode=ALL)
'''
    mark i as visited
    CLV = 0 
    SET = 0 
    while SET != 0 and CLV = 0:
        for v in SET:
            mark v as visited
            if L(v) == 1:
                CLV = Clv.append(v)
            else:
                for r in closest_neighbours(v):
                    if L(r) == 0:
                        SET = SET.append(r)
    SET_curr = SET_next
...
'''
    
    
    #closest_neighbours = assembly_graph.neighbors(i, mode=ALL) 
