#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import itertools
from operator import itemgetter
import sys
import numpy as np
import os

sampleid = sys.argv[1]
n_callers_filter_threshold = 2
threshold = 0.5

data = []

def loadme(filename):
    rows = []
    with open(filename) as f:
        for i in f:
            rows.append(i.strip().split("\t"))

    for i in range(len(rows)):
        rows[i][1] = int(rows[i][1])
        rows[i][2] = int(rows[i][2])

    return rows

for myfile in ["./" + sampleid + "_" + name + ".txt" for name in 
               ["cnvnator", "breakdancer", "lumpy","tardis", "pindel"]]:
    data.append(loadme(myfile))

chromos = []

for h,d in enumerate(data):
    chrom = {}
    start = None
    cur = None
    for k,r in enumerate(d):
        r[0] = int((r[0]).replace('X','23').replace('Y','24').replace('M','25'))
        if r[0] != cur:
            if start != None:
                chrom[cur] = [start,k]
            start = k
            cur = r[0]
    chrom[cur] = [start,k+1]
    chromos.append(chrom)

file_indexes = [i for i in range(len(data))]

file_combinations = [c for c in itertools.combinations(file_indexes,2)]

def match(r1, r2):
    chrom1 = r1[0]
    start1 = r1[1]
    end1 = r1[2]
    svlen1 = end1 - start1
    sv1 = r1[3]
    chrom2 = r2[0]
    start2 = r2[1]
    end2 = r2[2]
    svlen2 = end2 - start2
    sv2 = r2[3]
    if (chrom1 != chrom2): 
        return (0.0, chrom1, min(start1,start2), r1, r2, max(end1, end2), sv1, sv2)
    if start1 > end2 or start2 > end1:
        return (0.0, chrom1, min(start1,start2), r1, r2, max(end1, end2), sv1, sv2)
    if (sv1 == sv2) or (sv1 == 'MIXED') and (sv2 in ['DEL', 'DUP']) or (sv2 == 'MIXED') and (sv1 in ['DEL', 'DUP']):
        overlap = min(end1,end2) - max(start1,start2)    
        overlapmin = overlap/float(min(svlen1, svlen2))
        overlapmax = overlap/float(max(svlen1, svlen2))
        return (min(round(overlapmin,3), round(overlapmax,3)), chrom1, min(start1,start2), r1, r2, max(end1, end2), sv1, sv2)
    return (0.0, chrom1, min(start1,start2), r1, r2, max(end1, end2), sv1, sv2)

matches = []

for c in file_combinations:
    rows1 = data[c[0]]
    rows2 = data[c[1]]
    for i in range(len(rows1)):
        firstmatch = False
        r1 = rows1[i]
        if r1[0] not in chromos[c[1]]: continue
        start, end = chromos[c[1]][r1[0]] 

        for j in range(start, end):
            r2 = rows2[j]
            res = match(r1, r2)
            if res[0] >= threshold:
                matches.append(res[1:])
       
            if r2[0] > r1[0] or (r1[0] == r2[0] and r2[1] > r1[2]):
                break
                
matches_prefilter = matches
matches_sorted = sorted(matches_prefilter, key = itemgetter(0, 1))
matches_sorted = [set(map(tuple,m[2:4])) for m in matches_sorted]

matches_filtered = []
num_matches = len(matches_sorted)

for i in range(num_matches):
    m1 = matches_sorted[i]
    if i == num_matches - 1:
        matches_filtered.append(list(m1))
        break
    m2 = matches_sorted[i+1]
    if not m1.isdisjoint(m2):
        matches_sorted[i+1] = m1.union(m2)
    else:
        matches_filtered.append(list(m1))

matches_filtered_all = []
for i in matches_filtered:
    row_indexes = [j for j in range(len(i))]
    rows_combinations = [c for c in itertools.combinations(row_indexes,2)]
    all_local = []
    for r in rows_combinations:        
        rows1 = i[r[0]]
        rows2 = i[r[1]]
        res = match(rows1, rows2)
        all_local.append(res)
    matches_filtered_all.append(all_local)

below_threshold = 0
multi_software = 0

matches_good = []
matches_below = []
matches_multi = []

for i, row in enumerate(matches_filtered_all):
    ranges = [k[3][1:6] for k in row] + [k[4][1:6] for k in row]
    software = [k[3] for k in ranges]
    thresholds = [k[0] for k in row]

    good = True

    if len(set(ranges)) != len(set(software)):
        good = False
        multi_software +=1
        matches_multi.append(row)

    if any([m < threshold for m in thresholds]):
        good = False
        below_threshold += 1
        matches_below.append(row)

    if good:
        matches_good.append(row)

remove_indices = []

for index, row in enumerate(matches_good):
    software = []
    n_tuples = len(row)

    for n in range(n_tuples):
        software.append(row[n][3][4])
        software.append(row[n][4][4])

    n_software = len(set(software))
    if n_software < int(n_callers_filter_threshold):
        remove_indices.append(index)

matches_good = [i for j, i in enumerate(matches_good) if j not in remove_indices]


def magicfunction(mydict, mylist):
    for k in mylist:
        if k in mydict:
            return (mydict[k], k)

one_from_many = []

def most_frequent(List):
    return max(set(List), key = List.count)

for i, row in enumerate(matches_good):
    allrows = list(set([k[3] for k in row] + [k[4] for k in row]))
    posdict = {}
    genodict = {}
    chromosome = allrows[0][0]
    for r in allrows:
        if r[4] not in posdict: posdict[r[4]] = r[1:4]
        elif posdict[r[4]] != r[1:3]: print("unexpected pos")
        if r[4] not in genodict: genodict[r[4]] = r[5]
        elif genodict[r[4]] != r[3]: print("unexpected gt")

    if allrows[0][3] == "DEL":
        pos_winner, mypos = magicfunction(posdict,["lumpy","pindel","genomestrip","breakdancer","tardis","cnvnator"]) 
        gen_winner, mygen = magicfunction(genodict,["cnvnator","lumpy","genomestrip","tardis","breakdancer","pindel"])
        one = (chromosome, pos_winner, gen_winner, mypos, "_", mygen, r[6])
        one_from_many.append(one)

    if allrows[0][3] == "DUP":
        pos_winner, mypos = magicfunction(posdict,["pindel","lumpy","cnvnator","genomestrip","tardis","breakdancer"])
        gen_winner, mygen = magicfunction(genodict,["cnvnator","lumpy","genomestrip","tardis","breakdancer","pindel"])
        one = (chromosome, pos_winner, gen_winner, mypos, "_", mygen, r[6])
        one_from_many.append(one)
    
    # if allrows[0][3] == "INV":
    #     startS = [pos[1] for pos in allrows]
    #     endS = [pos[2] for pos in allrows]
    #     GTS = [pos[5] for pos in allrows]
    #     callers = '_'.join(list(set([pos[4] for pos in allrows])))
    #     start = round(np.median(startS))
    #     end = round(np.median(endS))
    #     GT = most_frequent(GTS)
    #     pos_winner = (start, end, 'INV')
    #     one = (chromosome, pos_winner, GT, callers, r[6])
    #     one_from_many.append(one)

    # if allrows[0][3] == "INS":
    #     startS = [pos[1] for pos in allrows]
    #     endS = [pos[2] for pos in allrows]
    #     GTS = [pos[5] for pos in allrows]
    #     callers = '_'.join(list(set([pos[4] for pos in allrows])))
    #     start = round(np.median(startS))
    #     end = round(np.median(endS))
    #     GT = most_frequent(GTS)
    #     pos_winner = (start, end, 'INS')
    #     one = (chromosome, pos_winner, GT, callers, r[6])
    #     one_from_many.append(one)

with open("./" + sampleid + "_callers_merged.txt", 'w') as o:
    for r in one_from_many: o.write("\t".join(map(str,r))+'\n')