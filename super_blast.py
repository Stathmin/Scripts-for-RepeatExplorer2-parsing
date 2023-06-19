import glob as gl
import pandas as pd
import re as re
import itertools as it
import xlsxwriter as xl
import textwrap as txtw
import numpy as np
import argparse
from blast import blaster
from blast import apply_leven as lev
import glob

from collections import defaultdict

###ОТРЕЗКИ

def intersec(a, b):
    #Находит у двух отрезков факт пересечения или 'почти' пересечения (5 нуклеотидов разницы)
    a1, a2, b1, b2 = int(a.split('-')[0]), int(a.split('-')[1]), int(b.split('-')[0]), int(b.split('-')[1])
    one = (min(a1, a2), max(a1, a2))
    two = (min(b1, b2), max(b1, b2))
    left = one if min(one) <= min(two) else two
    right = two if max(two) >= max(one) else one
    intersec_fact = min(right) - max(left) < 5 or max(right) - max(left) == 0
    return intersec_fact, left, right

def grouping(L):
    #проходит через все отрезки списка, попарно смотрит пересечения и их объединяет
    L = list(set(L))
    FullCycle=len(L)<2
    while not FullCycle:
        pairs = it.combinations(L, 2)
        changes = False
        for a, b in pairs:
            intersec_fact, left, right = intersec(a, b)
            if intersec_fact:
                changes = True
                c = str(min(left))+'-'+str(max(right))
                #print(L, left, right, a, b, (min(left), max(right)))
                if c == a and b in L:
                    L.remove(b)
                if c == b and a in L:
                    L.remove(a)
                else:
                    L.remove(a) if a in L else False
                    L.remove(b) if b in L else False
                    L.append(c) if c not in L else False
                if len(L)==1:
                    FullCycle = True
                break
        if not changes:
            FullCycle = True
    return('\n'.join(L))

    ##ПАРСИНГИ

def get_cyphers():
    with open('cypher.csv', 'r') as file:
        #Шифровки для образцов
        cypherdict = defaultdict(None)
        content = {x.split(',')[0]:x.split(',')[1:] for x in file.read().split('\n')}
        cypherdict.update(content)
        return cypherdict


def get_tareans(path, index):
    repeats_paths = gl.glob(f'{path}/*TAREAN*')
    names, seqs = [], []
    for repeat_path in repeats_paths:
        with open(repeat_path, 'r') as file:
            splited_fastas = [*filter(None, file.read().split('>'))]
            for fasta in splited_fastas:
                lines = [*filter(None, fasta.split('\n'))]
                names.append(index+'_'+lines[0])
                seqs.append('\n'.join(txtw.wrap("".join(lines[1:]), width=80)))

    cls = list(map(lambda x: re.findall('CL\d+',x)[0], names))
    cls_zeroes = list(map(lambda x: "CL{:04d}".format(int(re.findall('\d+',x)[0])), cls))
    pics = []
    for x in cls_zeroes:
        glob_search = gl.glob(f'{path}/seqclust/clustering/clusters/dir_{x}/graph_layout.png')
        if glob_search:
            pics.append(glob_search[0])
        else:
            pics.append('./No_image.png')
    tareandf = pd.DataFrame({'cluster':names, 'seq':seqs, 'pic_path':pics, 'Cluster':cls})
    return tareandf


def get_copies(path):
    cltab = gl.glob(f'{path}/CLUSTER_TABLE.csv')[0]
    copydf = pd.read_csv(cltab, skiprows=5, delimiter='\t')

    with open(cltab, 'r') as file:
        lines=file.readlines()
        total_size=int(lines[4].split('\t')[1].strip())

    copydf['Cluster']=copydf['Cluster'].apply(lambda x: 'CL'+str(x))
    copydf['size, %']=copydf['Size'].apply(lambda x: 100*x/total_size)
    return copydf


def get_comparative_data(index, path):
    df = pd.read_csv(f'{path}/COMPARATIVE_ANALYSIS_COUNTS.csv', skiprows=2, delimiter='\t', nrows=500)
    local_indexi = df.columns.drop(['cluster','supercluster']).to_list()
    for local_index in local_indexi:
        df[f'{local_index}, %'] = (100*df[local_index]/(df[local_indexi].sum(axis=1))).map('{:,.2f}'.format)
    df['merger'] = df['cluster'].apply(lambda x: index + '_CL' + str(x))
    return df.drop(['cluster', 'supercluster'], axis=1)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('index')

    args = parser.parse_args()

    indexi = [args.index]

    dbs = ['comp'] if len(re.findall('[A-Z]+[0-9]+', args.index))>1 else ['local']

    my_df = blaster(dbs=dbs, subjects=indexi, tasks=['megablast', 'dc-megablast', 'blastn'], verbose=False, num_threads=300)
    my_df.to_csv(f"{indexi[0]}_blast_result.csv")
