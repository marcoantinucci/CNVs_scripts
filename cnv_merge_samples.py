#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import itertools
from operator import itemgetter
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import os
import pickle
from sklearn import preprocessing

threshold = 0.5

data = []

def loadme(filename):
    rows = []
    with open(filename) as f:
        for i in f:
            rows.append(tuple(i.strip().split("\t")))
    return rows

listsamples = []
with open("./samples_list.txt") as s:
    for i in s:
        listsamples.append(i.strip())


for myfile in ["./" + sampleid + "_callers_merged.txt" for sampleid in listsamples]:
    data.append(loadme(myfile))

field_3_encoder = preprocessing.LabelEncoder()
field_3_encoder.fit(['DEL','DUP'])
field_4_encoder = preprocessing.LabelEncoder()
field_4_encoder.fit(['0/0', '0/1', '1/1'])
field_5_encoder = preprocessing.LabelEncoder()
sw = ['cnvnator', 'breakdancer', 'pindel', 'tardis', 'lumpy', 'genomestrip']
field_5_encoder.fit(['_'.join(s) for s in itertools.product(sw,sw)])
field_6_encoder = preprocessing.LabelEncoder()
field_6_encoder.fit(listsamples)

dt=[('0','int64'),('1','int64'),('2','int64'),('3','<U3'),('4','<U3'),('5','<U32'),('6','<U32')]   
dt_=[('0','uint8'),('1','int32'),('2','int32'),('3','uint8'),('4','uint8'),('5','uint8'),('6','uint16'),('7','int32')]

standard = {'0':'uint8','1':'uint32','2':'uint32','3':'uint8','4':'uint8','5':'uint8','6':'uint16'}
standard_p = {'0':'uint8','1':'uint32','2':'uint32','3':'uint8','4':'uint8','5':'uint8','6':'uint16','7':'int32'}

data_ = []
for i, d in enumerate(data):
    xarr = np.array(d,dtype=dt) 
    xarr_ = np.empty((len(xarr),), dtype=dt_)
    xarr_['0'] = xarr['0'].astype('uint8')
    xarr_['1'] = xarr['1'].astype('int32')
    xarr_['2'] = xarr['2'].astype('int32')
    xarr_['3'] = field_3_encoder.transform(xarr['3']).astype('uint8')
    xarr_['4'] = field_4_encoder.transform(xarr['4']).astype('uint8')
    xarr_['5'] = field_5_encoder.transform(xarr['5']).astype('uint8')
    xarr_['6'] = field_6_encoder.transform(xarr['6']).astype('uint16')
    xarr_['7'] = ((xarr['1'] + xarr['2']) / np.int32(2)).astype('int32')

    df=pd.DataFrame(xarr_)
    df = df[df['4'] != 0]
    data_.append(df.set_index(['3'], drop=False).sort_values('7'))

del data
data = data_

file_indexes = range(len(data))

out = []
times = 0
totj = len(file_indexes)     

allr = pd.concat(data, copy=True)
allr.sort_values('7', inplace=True)
enum = field_6_encoder.transform(listsamples).astype('uint16')

for i in file_indexes:
    times += 1    
    rows1 = data[i]
    rows2 = allr[allr['6'] != enum[i]]
    for sv in [np.uint8(0), np.uint8(1)]:
        m = pd.merge_asof(rows1[rows1.index==sv], rows2[rows2.index==sv], on='7', by='0', tolerance=np.int32(2**16), direction='nearest') 
        m.dropna(inplace=True)

        m = m[(m['1_x'] < m['2_y']) & (m['1_y'] < m['2_x'])]
        svlen1 = m['2_x'] - m['1_x']
        svlen2 = m['2_y'] - m['1_y']
        m['overlap'] = (m[['2_x','2_y']].min(axis=1) - m[['1_x','1_y']].max(axis=1)) / pd.concat([svlen1, svlen2], axis=1).max(axis=1)
        m = m[m['overlap']>=threshold]

        if len(m):
            cols = ['0','1','2','3','4','5','6']
            mx = m[['0','1_x','2_x','3_x','4_x','5_x','6_x']]
            mx.columns = cols
            my = m[['0','1_y','2_y','3_y','4_y','5_y','6_y']]
            my.columns = cols
            my = my.astype(standard)
            m = pd.concat((mx, my), copy=False)
            m.drop_duplicates(inplace=True)
            out.append(m)

n = pd.concat(out, copy=False)
n.drop_duplicates(inplace = True)

n['7'] = np.floor_divide(n['2'] + n['1'], np.uint32(2)).astype('uint32')

n = n.set_index(['3'], drop=False).sort_values(['0','7'])

n_ = n
ofin = []
maxindex = np.uint32(0)

for sv in [np.uint8(0), np.uint8(1)]:
    n = n_[n_.index==sv]
    ns = n.shift(-1)
    ns.dropna(inplace=True)
    ns = ns.astype(standard_p)
    ns.columns = [c + '_s' for c in n.columns]
    sh = pd.concat((n[:-1],ns), axis=1)
    sh['overlap'] = (sh[['2','2_s']].min(axis=1) - sh[['1','1_s']].max(axis=1)) / pd.concat([sh['2'] - sh['1'], sh['2_s'] - sh['1_s']], axis=1).max(axis=1)
    sh['overlap'] = (sh['0'] == sh['0_s']) * sh['overlap']
    sh['overlap'] = (sh['overlap']<threshold).astype('uint8')
    sh['k'] = (sh['overlap']>threshold).cumsum().astype('uint32') + maxindex
    sh = sh[sh.overlap==0]
    cols = ['0','1','2','3','4','5','6','k']
    mx = sh[['0','1','2','3','4','5','6', 'k']]
    mx.columns = cols
    my = sh[['0_s','1_s','2_s','3_s','4_s','5_s','6_s','k']]
    my.columns = cols
    fin = pd.concat((mx, my), copy=False)
    fin.drop_duplicates(inplace=True)
    fin.sort_values(['0','1'], inplace=True)
    maxindex += fin.k.max() + np.int32(1)
    ofin.append(fin)

fin = pd.concat(ofin, copy=False)

singletones_control = sys.argv[1]
if singletones_control == "yessingletones":
    fin.reset_index(drop=True, inplace=True)
    mc = ['0','1','2','3','4','5','6']
    left = pd.merge(allr.reset_index(drop=True)[mc], fin, on=mc, how='left')
    singletons = left[left.k.isnull()][mc]
    fin['3'] = field_3_encoder.inverse_transform(fin['3'])
    fin['4'] = field_4_encoder.inverse_transform(fin['4'])
    fin['5'] = field_5_encoder.inverse_transform(fin['5'])
    fin['6'] = field_6_encoder.inverse_transform(fin['6'])
    singletons['3'] = field_3_encoder.inverse_transform(singletons['3'])
    singletons['4'] = field_4_encoder.inverse_transform(singletons['4'])
    singletons['5'] = field_5_encoder.inverse_transform(singletons['5'])
    singletons['6'] = field_6_encoder.inverse_transform(singletons['6'])
    singletons.sort_values(['0','1'], inplace=True)
    singletons['k'] = -1
    fin = pd.concat([fin, singletons], ignore_index=True)
    fin.sort_values(['0','1'], inplace=True)
    fin.reset_index(drop = True, inplace=True)
    
    firstsingletone = False
    for index, row in fin.iterrows():
        k = row[7]
        if index == len(fin) - 1:
            if k == -1:
                last_k += 1
                fin.at[index, 'k'] = last_k
            else:
                prev_k = fin.at[index - 1, 'k']
                if k == prev_k:
                    fin.at[index, 'k'] = last_k
                else:
                    last_k += 1
                    fin.at[index, 'k'] = last_k                
            break
        next_k = fin.at[index + 1, 'k']
        if index == 0 and k == -1:
            fin.at[index, 'k'] = 1
            last_k = 2
            firstsingletone = True
        if index == 0 and firstsingletone == False:
            last_k = 1
        if index > 0:
            if k == -1:
                last_k += 1
                fin.at[index, 'k'] = last_k
                if next_k != -1:
                    last_k += 1
            else:
                if k == next_k:
                    fin.at[index, 'k'] = last_k
                elif next_k != -1:
                    fin.at[index, 'k'] = last_k
                    last_k += 1
                elif next_k == -1:
                    fin.at[index, 'k'] = last_k


if singletones_control == "nosingletones":
    fin['3'] = field_3_encoder.inverse_transform(fin['3'])
    fin['4'] = field_4_encoder.inverse_transform(fin['4'])
    fin['5'] = field_5_encoder.inverse_transform(fin['5'])
    fin['6'] = field_6_encoder.inverse_transform(fin['6'])


dataset_groups = pd.read_csv("./samples_study.txt", sep = "\t")

fin2 = fin.merge(dataset_groups[['sampleID', 'study']], left_on='6', right_on='sampleID').drop('sampleID', axis='columns')

if singletones_control == "yessingletones":
    df = pd.DataFrame(fin2)
    df.to_csv("./samples_merged_flat.txt", index = False, header = False, sep="\t")
if singletones_control == "nosingletones":
    df = pd.DataFrame(fin2)
    df.to_csv("./samples_merged_NOSINGLE_flat.txt", index = False, header = False, sep="\t")

def formatter(common_variants):
    starts = []
    ends = []
    samplesstring = []
    samples = []
    result = ()
    for index in range(len(common_variants)):
        singletone = False
        shared_SV = common_variants[index][7]
        starts.append(int(common_variants[index][1]))
        ends.append(int(common_variants[index][2]))
        samples.append(common_variants[index][6])
        samplesstring.append('_'.join(common_variants[index][0:7]))
        if index == len(common_variants) -1:
            if shared_SV in k_singletones:
                singletone = True
                startline = common_variants[index][0:4]
                sample = common_variants[index][6]
                samplestring = '_'.join(common_variants[index][0:7])
                samplesstring = []
                samples = []
                starts = []
                ends = []
            else:
                minstart = min(starts)
                maxend = max(ends)
                startline = common_variants[index][0:4]
                startline[1] = minstart
                startline[2] = maxend
            for name in listsamples:
                if not singletone:
                    if name in samples:
                        samples2 = [x.split("_")[7] for x in samplesstring]
                        selectindex = samples.index(name)
                        selectstring = [s for s in samplesstring if name in s]
                        startline.append(samplesstring[selectindex])
                    else:
                        startline.append('0/0')
                if singletone:
                    if name != sample:
                        startline.append('0/0')
                    else:
                        startline.append(samplestring)
            break    
        if shared_SV == common_variants[index + 1][7] and shared_SV not in k_singletones:
            continue
        else:
            if shared_SV in k_singletones:
                singletone = True
                startline = common_variants[index][0:4]
                sample = common_variants[index][6]
                samplestring = '_'.join(common_variants[index][0:7])
                samplesstring = []
                samples = []
                starts = []
                ends = []
            else:
                minstart = min(starts)
                maxend = max(ends)
                startline = common_variants[index][0:4]
                startline[1] = minstart
                startline[2] = maxend

        for name in listsamples:
            if not singletone:
                if name in samples:
                    samples2 = [x.split("_")[7] for x in samplesstring]
                    selectindex = samples.index(name)
                    selectstring = [s for s in samplesstring if name in s]
                    startline.append(samplesstring[selectindex])
                else:
                    startline.append('0/0')
            if singletone:
                if name != sample:
                    startline.append('0/0')
                else:
                    startline.append(samplestring)
        samplesstring = []
        samples = []
        starts = []
        ends = []
        result += (startline,)

    df = pd.DataFrame(result)
    df.to_csv("./samples_merged.txt", index = False, header = False, sep = "\t")


input_variants = []
with open("./samples_merged_flat.txt") as f:
    for k in f:
        input_variants.append(k.strip().split("\t"))

formatter(input_variants)
