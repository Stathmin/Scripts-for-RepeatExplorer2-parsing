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

###sections

def intersec(a, b):
    #finds intersection or near-intersection (up-to 5 bp)
    a1, a2, b1, b2 = int(a.split('-')[0]), int(a.split('-')[1]), int(b.split('-')[0]), int(b.split('-')[1])
    one = (min(a1, a2), max(a1, a2))
    two = (min(b1, b2), max(b1, b2))
    left = one if min(one) <= min(two) else two
    right = two if max(two) >= max(one) else one
    intersec_fact = min(right) - max(left) < 5 or max(right) - max(left) == 0
    return intersec_fact, left, right

def grouping(L):
    #iterates over all the sections, checks for intersections, unites if necessary
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

def get_cyphers():
    with open('cypher.csv', 'r') as file:
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
    local_indexes = df.columns.drop(['cluster','supercluster']).to_list()
    for local_index in local_indexes:
        df[f'{local_index}, %'] = (100*df[local_index]/(df[local_indexes].sum(axis=1))).map('{:,.2f}'.format)
    df['merger'] = df['cluster'].apply(lambda x: index + '_CL' + str(x))
    return df.drop(['cluster', 'supercluster'], axis=1)

def get_columns(x):
    try:
        return(list(x.items())[-1][-1].columns.values.tolist())
    except:
        return([])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('index')
    parser.add_argument('local_or_comparative')

    args = parser.parse_args()

    indexes = [args.index]

    all = set(['KK1', 'KK3', 'KK6', 'KA25', 'KA27']) #indexes list, from each to deduct which are not align with the repeat
    tasks=['megablast', 'dc-megablast', 'blastn']

    dbs = [args.local_or_comparative]
    paths = [f'./{x}' for x in indexes]

    my_df = blaster(dbs=dbs, subjects=indexes, tasks=tasks, verbose=False, num_threads=300)
    ncbi_df = my_df[my_df['db'].apply(lambda x: 'ncbi' in x)]
    local_df = my_df[my_df['db'].apply(lambda x: ('local' in x or 'comp' in x))]
    ref_df = my_df[my_df['db'].apply(lambda x: x=='ref')]
    ind_path_dict = dict(zip(indexes, paths))

    dataframe_list = list()

    with open('./ncbi_naming.csv','r') as file:
        ncbi_naming = pd.read_csv(file).drop('Unnamed: 0', axis=1)

    for index, path in ind_path_dict.items():
        print(index, path)
        df_tar = get_tareans(path, index)
        df_cop = get_copies(path)
        df = pd.merge(df_tar, df_cop, left_on='Cluster', right_on='Cluster').drop(['Cluster', 'Final_annotation'], axis=1)
        if dbs == ['comp']:
            df_comp = get_comparative_data(index, path)
            df['merger'] = df['cluster'].apply(lambda x: "_".join(x.split('_')[0:2]))
            df = pd.merge(df, df_comp, left_on='merger', right_on='merger').drop(['merger'], axis=1)
            df['comparative']=True
        dataframe_list.append(df)

    annotation_db = pd.concat(dataframe_list)
    if 'comparative' in annotation_db.columns.values.tolist():
        annotation_db['comparative'] = annotation_db['comparative'].fillna(False)
    annotation_db = annotation_db[annotation_db['TAREAN_annotation'] != 'rDNA']
    annotation_db = annotation_db.rename({'size, %':'proportion, %', 'Automatic_annotation':'TAREAN annotation'}, axis=1)

    ncbi_df = my_df[my_df['db'].apply(lambda x: 'ncbi' in x)]
    local_df = my_df[my_df['db'].apply(lambda x: ('local' in x) or ('comp' in x))]
    ref_df = my_df[my_df['db'].apply(lambda x: x=='ref')]

    tarean_dict = {}
    for type_gr, group in annotation_db.groupby('TAREAN_annotation'):
        group['proportion, %'] = group['proportion, %'].apply(lambda x:'{0:.3f}'.format(x))
        percent_columns = [col for col in group.columns if '%' in col]
        group = group[['cluster', 'seq', 'pic_path', 'proportion, %', 'TAREAN annotation']+percent_columns]
        tarean_dict.update({type_gr: group})

    ncbi_dict = {}
    for qseqid, group in ncbi_df.groupby('qseqid'):
        group = group.sort_values(['task', 'fitting', 'sseqid'])[group['sum_coverage']>0.7]
        group['db'] = 'ncbi'
        if len(group) > 0:
            group = lev(group).head()
        group = group.merge(ncbi_naming, left_on='sseqid', right_on='short').drop('sseqid', axis=1)[['full', 'task', 'sum_coverage']]
        group['sum_coverage'] = group['sum_coverage'].apply(lambda x:'{0:.3f}'.format(x))
        group = group.rename({'full':'homology to repeat-related sequences in Nucleotide DB','sum_coverage':'share of cluster consensus covered'}, axis=1)
        group = group[['homology to repeat-related sequences in Nucleotide DB', 'task', 'share of cluster consensus covered']]
        ncbi_dict.update({qseqid:group})

    cypherdict = get_cyphers()
    lmb = lambda x: x.split(' ')[0][0]+'.'+x.split(' ')[1]
    SHDICT= defaultdict(str)
    SHDICT.update(dict([(x,lmb(y[0])) for x, y in cypherdict.items()]))

    local_dict = {}
    for qseqid, group in local_df.groupby('qseqid'):
        group = group.sort_values(['task', 'fitting', 'sseqid'])[group['sum_coverage']>0.7]
        group['db'] = 'local'
        if len(group) > 0:
            group['sseqid'] = group['sseqid'].apply(lambda x: SHDICT[x.split('_')[0].split('(')[0]]+' '+x)
            group['sum_coverage'] = group['sum_coverage'].apply(lambda x:'{0:.3f}'.format(x))
        group = group.rename({'sseqid':'homology to other sequences in this work','sum_coverage':'share of cluster consensus covered'}, axis=1)
        group = group[['homology to other sequences in this work', 'task', 'share of cluster consensus covered']]
        found = set(group['homology to other sequences in this work'].apply(lambda x: re.findall(' .*_CL', x)[0].replace('_CL','').replace(' ','')).tolist())
        not_found = all-found-set([index, 'KK1KK2KK3KK4'])
        if len(not_found) > 0:
            group.loc[len(group)] = ["*"+",".join(sorted(list(not_found))),'blastn','NOT FOUND']
        local_dict.update({qseqid:group})

    ref_dict = {}
    for qseqid, group in ref_df.groupby('qseqid'):
        group = group.sort_values(['task', 'sseqid'])
        group = group.rename({'sseqid':'Homology to oligos and repeats from other works','sum_coverage':'share of cluster consensus covered'}, axis=1)
        group = group[['Homology to oligos and repeats from other works','pident', 'length', 'evalue']]
        ref_dict.update({qseqid:group})

    workbook = xl.Workbook(f'./solo_annotation/{indexes[0]}.xlsx')
    worksheet = workbook.add_worksheet('main_sheet')
    sidecols=2

    columnsets = (get_columns(tarean_dict),get_columns(local_dict),get_columns(ncbi_dict),get_columns(ref_dict))
    sets = ['annotation', 'local', 'ncbi', 'ref']
    coords = {i:sum([len(x) for x in columnsets[0:num]])+sidecols for num, i in enumerate(sets)}
    columnnames = [x for l in columnsets for x in l ]

    for num, i in enumerate(columnnames):
        worksheet.write(0,num+sidecols,i)


    merge_format = workbook.add_format({
    'bold': 1,
    'border': 1,
    'align': 'center',
    'valign': 'vcenter'})

    merge_format_seq = workbook.add_format({
    'bold': 1,
    'border': 1,
    'align': 'left',
    'valign': 'vcenter'})


    row = 1
    type_order = ['Putative satellites (high confidence)', 'Putative satellites (low confidence)', 'Putative LTR elements']
    for key in type_order:
        worksheet.merge_range(row, min(coords.values()), row, max(coords.values())+3, key, merge_format)
        row+=1
        if key in tarean_dict.keys():
            annotation = tarean_dict[key]
        else:
            continue
        for _, an_row in annotation.iterrows():
            cluster = an_row['cluster']
            annotation_chunk = an_row.to_frame().transpose()
            local_chunk = local_dict[cluster] if cluster in local_dict.keys() else list()
            ncbi_chunk = ncbi_dict[cluster] if cluster in ncbi_dict.keys() else list()
            ref_chunk = ref_dict[cluster] if cluster in ref_dict.keys() else list()
            chunks = [annotation_chunk, local_chunk, ncbi_chunk, ref_chunk]
            set_chunks = dict(zip(sets, chunks))

            for position, chunk in set_chunks.items():
                col = coords[position]
                if len(chunk) > 0:
                    for rel_row, data in chunk.reset_index(drop='index').iterrows():
                        for rel_col, item in enumerate(data):
                            if position != 'annotation':
                                if 'ImportantIndex_' in str(item): #example of style marking for readability
                                    worksheet.write(row+rel_row, col+rel_col, item, workbook.add_format({'bold': True, 'font_color': 'red'}))
                                else:
                                    worksheet.write(row+rel_row, col+rel_col, item)
                            else:
                                if rel_col in [2]: #to anchor graph layout pictures
                                    path = item
                                    item = ''
                                if row+max(list(map(len, chunks)))-1 != row:
                                    worksheet.merge_range(row, col+rel_col, row+max(list(map(len, chunks)))-1, col+rel_col, item, merge_format)
                                else:
                                    worksheet.write(row, col+rel_col, item, merge_format)
                                if rel_col in [2]:
                                    worksheet.insert_image(row, col+rel_col, path, {'x_scale': 0.15, 'y_scale': 0.15, 'object_position': 1})
            row = row + max(list(map(len, chunks)))

    worksheet.write(1, 0, 'Paired reads')
    worksheet.write(2, 0, '% reads in clusters')
    with open(gl.glob(f'{indexes[0]}/index.html')[0], 'r') as file:
        text = file.read()
        pair_number_context = re.findall("<p class='character'>Number of analyzed reads: [0-9]+</p>", text)[0]
        pair_number_an = int(re.findall('[0-9]+', pair_number_context)[0])
        percent_number_context = re.findall("<p class='character'>Proportion of reads in top clusters : [0-9]+ %</p>", text)[0]
        percent_number = int(re.findall('[0-9]+', percent_number_context)[0])
    worksheet.write(1, 1, pair_number_an)
    worksheet.write(2, 1, percent_number)
    worksheet.merge_range(0,0,0,1, SHDICT[an_row['cluster'].split('_')[0].split('(')[0]], merge_format)

    worksheet.set_column(2,2, 30)
    worksheet.set_column(4,4, 40)
    worksheet.set_column(7,7, 40)
    worksheet.set_column(8,8, 40)
    worksheet.set_column(9,9, 40)
    worksheet.set_column(12,12, 40)
    workbook.close()