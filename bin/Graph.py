#!/usr/bin/env python3

import glob
import pandas as pd
import copy
import numpy as np
import warnings
import os
import sys

filter_threshold = 0.6


# return 1 if x is '+', -1 if x is '-', 0 otherwise
def sign(x):
    if x == '+':
        return 1
    if x == '-':
        return -1
    else:
        return 0


# read in the graph file and return two dataframes: segments and lines
# segments have five columns ['sign': 'S', 'contig_ID': int,'sequence':string, 'length':int,'depth': float]
# lines have five columns ['sign': 'L','start':int,'start_type': string, '+' or '-', 'end': int, 'end_type':string, '+' or '-']
# It also reads empty files and return empty dataframes
def readUniGraph(file_name):
    f = open(file_name)
    segments = pd.DataFrame(columns=['sign', 'contig_ID', 'sequence', 'length', 'depth'])
    lines = pd.DataFrame(columns=['sign', 'start', 'start_type', 'end', 'end_type'])
    for line in f.readlines():
        line = line.split()
        if line[0] == 'S':  # read segments
            s = copy.deepcopy(line)
            s[3] = line[3].split(":")[-1]
            s[4] = line[4].split(":")[-1]
            s_series = pd.Series(s, index=segments.columns)
            segments = segments.append(s_series, ignore_index=True)
        elif line[0] == 'L':  # read lines
            s = line[:5]
            s_series = pd.Series(s, index=lines.columns)
            lines = lines.append(s_series, ignore_index=True)

    # convert datatype
    segments.length = segments.length.astype(int)
    segments.contig_ID = segments.contig_ID.astype(int)
    segments.depth = segments.depth.astype(float)
    lines.start = lines.start.astype(int)
    lines.end = lines.end.astype(int)
    segments.index = segments.contig_ID.values
    # print(lines)
    lines['signed_start'] = lines.start_type.apply(sign) * lines.start
    lines['signed_end'] = lines.end_type.apply(sign) * lines.end
    return segments, lines


# delete a segment
def remove_segment(S, L, s):
    # remove the segment from segments table
    segments = S[S.contig_ID != s]
    # also remove the lines connected to those segments.
    lines = L[L.start != s]
    lines = lines[lines.end != s]
    return segments, lines


# remove non-connected small contigs
def filter_contigs(S, L):
    if not len(S):
        return S, L
    connected = L.start.unique().tolist() + L.end.unique().tolist()
    segments = pd.DataFrame()
    left = pd.DataFrame()
    for s in S.index:
        if s in connected:
            segments = segments.append(S.loc[s, :])
        if s not in connected:
            left = left.append(S.loc[s, :])
    if len(left):
        left['weight'] = left.length * left.depth
        left.sort_values(by=['weight'], ascending=False, inplace=True)
        segments = segments.append(left.head(3))
    segments.sort_index(inplace=True)
    return segments[S.columns], L


# remove ends with less coverage than the threshold
def filter_ends(S, L, threshold=filter_threshold):
    segments = copy.deepcopy(S)
    lines = copy.deepcopy(L)

    start = lines.groupby('start')['signed_start'].unique()
    end = lines.groupby('end')['signed_end'].unique()
    end = (-end)

    for i in end.index:
        if i in start:
            start[i] = np.array(start[i].tolist() + end[i].tolist())
        else:
            start[i] = end[i]

    for s in start.index:
        l = pd.Series(start[s]).nunique()
        if l == 1 and segments.depth[s] < threshold:
            segments, lines = remove_segment(segments, lines, s)

    return segments, lines


# read in two dataframe S (segments) and L(lines)
# read in the threshold of sequencing depth
# return two filtered dataframe, s(segments) and l(lines)
def second_filter(S, L, threshold=filter_threshold):
    if not len(S):
        return S, L
    segments = copy.deepcopy(S)
    lines = copy.deepcopy(L)
    # 1. filter small, non-connected contigs
    # segments, lines = filter_contigs(S, L)
    # 2. filter out end segments with less than average coverage
    average_coverage = sum(segments.depth * segments.length) / sum(segments.length)
    segments.depth = segments.depth / average_coverage
    segments, lines = filter_ends(segments, lines, threshold)

    # 3. filter coverage
    # repeat until the average coverage of all segments are above threshold
    average_coverage = sum(segments.depth * segments.length) / sum(segments.length)
    segments.depth = segments.depth / average_coverage
    segments.depth = segments.apply(lambda x: x.depth if (x.length > 10 or x.depth > threshold) else 1, axis=1)
    while sum(segments.depth < threshold * 0.8):
        segments = segments[segments.depth >= threshold * 0.8]  # modify with depth
        nodes = segments.contig_ID.unique()
        if len(lines):
            lines = lines[lines.start.apply(lambda x: True if x in nodes else False)]
            lines = lines[lines.end.apply(lambda x: True if x in nodes else False)]
        average_coverage = sum(segments.depth * segments.length) / sum(segments.length)
        segments.depth = segments.depth / average_coverage
        segments.depth = segments.apply(lambda x: x.depth if (x.length > 10 or x.depth > threshold) else 1, axis=1)

    # if no line left, return all the remaining segments and an empty DataFrame of lines
    if len(lines) == 0:
        lines = pd.DataFrame(columns=L.columns)
        return segments, lines

    # 4. filter self loops. Add or remove self-loop based on the sequencing depth.
    loops = lines[lines.start == lines.end]
    count_1 = abs(lines.start).value_counts()
    count_2 = abs(lines.end).value_counts()
    for i in loops.start:
        x = (count_1[i] + count_2[i]) / 2
        depth = round(segments.loc[i, 'depth'] - x)
        if depth >= 1:
            if depth > 5:  # control the largest depth to guarantee the running speed of largest().
                depth = 5
                warnings.warn("complex repeat region detected")  # raise a warning for heavily repeated regions.
            new_l = loops[loops.start == i]
            depth = int(depth)
            lines = lines.append([new_l] * depth, ignore_index=True)
        elif depth < -1:
            lines = lines[(lines.start != i) | (lines.end != i)]

    return segments, lines


# turn the dataframe(lines) to four lists in a fixed order s(start), e(end), rs(reverse_start), re(reverse_end)
# for each index(i), s[i]-->e[i] and rs[i]-->re[i] are the two directions of a line.
def Graph_to_list(L):
    lines = copy.deepcopy(L)
    lines['re_start'] = -lines.signed_end  # - end is reverse start
    lines['re_end'] = -lines.signed_start  # - start is the reverse end
    s = lines.signed_start.to_list()
    e = lines.signed_end.to_list()
    rs = lines.re_start.to_list()
    re = lines.re_end.tolist()
    return s, e, rs, re


# Return the reverse compliment of a DNA sequence
def reverse_compliment(s):
    y = ''
    for x in s:
        if (x == 'A') or (x == 'a'):
            y += 'T'
        if (x == 'T') or (x == 't'):
            y += 'A'
        if (x == 'C') or (x == 'c'):
            y += 'G'
        if (x == 'G') or (x == 'g'):
            y += 'C'
    return y[::-1]


# Read in the segments, largest list (L), and size of the largest list returned from find()
# Return two dataframes (s_list_to_graph, l_list_to_graph) containing the information of the best assembly
def write_result(S, Lists, filename):
    name = filename.split('.')[0]
    gfa_name = name + '.final.gfa'
    fasta_name = name + '.final.fasta'
    info_name = name + '.info.csv'

    if not Lists:  # empty output from last step.
        open(info_name, 'a').close()
        open(fasta_name, 'a').close()
        open(gfa_name, 'a').close()
        return 0

    if len(Lists) == 1:  # single segment without any connection
        start = abs(Lists[0])
        s = S.sequence[start]
        length = S.length[start]
        circle = False
        l = ''
        ambiguity = 1 - (s.count('A') + s.count('T') + s.count('C') + s.count('G')) / len(s)
    elif len(Lists) >= 2:
        lists = copy.deepcopy(Lists)
        s = ''
        l = ''
        circle = False
        start = lists[0]
        if lists[0] == lists[-1]:  # circularized assembly
            lists.pop(-1)
            circle = True
            l = 'L\t%d\t+\t%d\t+\t0M\n' % (abs(start), abs(start))
        for i in lists:  # adding each segments
            if i > 0:  # forward connection
                x = S.sequence[i]
            else:  # reverse connection
                x = S.sequence[-i]
                x = reverse_compliment(x)
            s += x
        length = len(s)
        ambiguity = 1 - (s.count('A') + s.count('T') + s.count('C') + s.count('G')) / len(s)

    # save info
    info = pd.DataFrame(columns=['ID', 'Length', 'Circle', 'Ambiguity'], data=[[name, length, circle, ambiguity]])
    info.to_csv(info_name, index=False)

    # save gfa
    g = 'S\t%d\t%s\tLN:i:%d\tdp:f:1.0' % (abs(start), s, length) + '\n' + l
    with open(gfa_name, 'w+') as text_file:
        text_file.write(g)

    # save fasta
    f = '>1\tlength=%d\tdepth=1.0x\tcircle=%r\n%s\n' % (length, circle, s)
    with open(fasta_name, 'w+') as text_file:
        text_file.write(f)
    return 1


# find all the index in a list(seq)
# return a list of the index
def search_index(x, seq):
    # read in two dataframes: segments and lines, filter by depth, return lists of lines
    ind = []
    for i in range(len(seq)):
        if x == seq[i]:
            ind.append(i)
    return ind


# Read in the information of segments, lists(from Graph_to_list), two numbers (start, end)
# find the longest connected segments that start with 'start'
# if the resulting connection is circularized, adding additional value to it.
# return a list of segments for the longest connection and its size
def largest(segments, lists, start, end):
    start_list, end_list, re_start_list, re_end_list = copy.deepcopy(lists)
    if len(start_list) == 0:  # empty lists
        return [], 0
    else:
        x1 = search_index(start, start_list)
        x2 = search_index(start, re_start_list)
        s1 = 0
        s2 = 0
        l1 = []
        l2 = []
        largest_size = 0
        largest_list = []
        # search through all index
        for x in x1:  # forward direction
            if end_list[x] == end:  # circularized and ended.
                l1 = [end_list[x]]
                s1 = 1000 + segments.length[abs(end)] / 1000
            # search recursively
            new_list = [start_list[:x] + start_list[x + 1:], end_list[:x] + end_list[x + 1:],
                        re_start_list[:x] + re_start_list[x + 1:], re_end_list[:x] + re_end_list[x + 1:]]
            l, s = largest(segments, new_list, end_list[x], end)
            size = s + segments.length[abs(end_list[x])] / 1000
            if size >= largest_size:
                largest_size = size
                largest_list = [end_list[x]] + l
            if s1 > largest_size:
                largest_size = s1
                largest_list = l1
        for x in x2:  # reverse direction
            if re_end_list[x] == end:
                l2 = [re_end_list[x]]
                s2 = 1000 + segments.length[abs(end)] / 1000
            # search recursively
            new_list = [start_list[:x] + start_list[x + 1:], end_list[:x] + end_list[x + 1:],
                        re_start_list[:x] + re_start_list[x + 1:], re_end_list[:x] + re_end_list[x + 1:]]
            l, s = largest(segments, new_list, re_end_list[x], end)
            size = s + segments.length[abs(re_end_list[x])] / 1000
            if size >= largest_size:
                largest_size = size
                largest_list = [re_end_list[x]] + l
            if s2 > largest_size:
                largest_size = s2
                largest_list = l2
    # print(largest_size, start, largest_list)
    return largest_list, largest_size


# loop for start point, find the longest connected map of a graph (in segments and lines)
# return a list of segments for the longest connection and its size
def find(S, L):
    segments, lines = second_filter(S, L)
    nodes = list(segments.contig_ID.unique())  # get a list of all the segments
    nodes.sort()
    largest_list = []
    largest_size = 0
    lists = Graph_to_list(lines)
    # print(nodes)
    # print(lines)
    while len(nodes):  # loop all the nodes
        start = nodes.pop(0)
        l1, s1 = largest(segments, lists, start, start)
        l2, s2 = largest(segments, lists, -start, -start)
        # calculate the length and add weights for circle
        s1 = (s1 % 1000 * 1000 + segments.length[start]) * (1 + s1 // 1000)
        s2 = (s2 % 1000 * 1000 + segments.length[start]) * (1 + s2 // 1000)

        if s1 > largest_size:
            largest_size = s1
            largest_list = [start] + l1
            if l1 and start == l1[-1]:
                nodes = list(set(nodes) - set([abs(x) for x in l1]))  # move visited nodes
        if s2 > largest_size:
            largest_size = s2
            largest_list = [-start] + l2
            if l2 and start == l2[-1]:
                nodes = list(set(nodes) - set([abs(x) for x in l2]))
    #      print(start, largest_list)
    return largest_list, largest_size


paths = [path for path in sorted(os.listdir('.')) if path.endswith('.gfa')]
# paths = ['/Users/liwang/Documents/plasmid/Graph/SQW-0072_SO10875_LR/SQW-0072_SO10875_B01.gfa']
for f in paths:
    name = f.split('/')[-1].split('.')[-2]
    # 1. read the file
    original_segments, original_lines = readUniGraph(f)
    segments, lines = filter_contigs(original_segments, original_lines)
    nodes_1, size_1 = find(segments, lines)
    write_result(original_segments, nodes_1, name)
    for node in nodes_1:
        segments, lines = remove_segment(segments, lines, node)
    nodes_2, size_2 = find(segments, lines)
#    if size_2:
#        write_result(original_segments, nodes_2, name + '_2')
