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

    paths = [f'./{x}' for x in indexi]

    my_df = blaster(dbs=dbs, subjects=indexi, tasks=['megablast', 'dc-megablast', 'blastn'], verbose=False, num_threads=300)
    ncbi_df = my_df[my_df['db'].apply(lambda x: 'ncbi' in x)]
    local_df = my_df[my_df['db'].apply(lambda x: ('local' in x or 'comp' in x))]
    ref_df = my_df[my_df['db'].apply(lambda x: x=='ref')]

    ind_path_dict = dict(zip(indexi, paths))

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
            print(df_comp.columns)
            df['merger'] = df['cluster'].apply(lambda x: "_".join(x.split('_')[0:2]))
            print(df.columns)
            df = pd.merge(df, df_comp, left_on='merger', right_on='merger').drop(['merger'], axis=1)
            df['comparative']=True
        dataframe_list.append(df)

    annotation_db = pd.concat(dataframe_list)
    print(annotation_db)
    if 'comparative' in annotation_db.columns.values.tolist():
        annotation_db['comparative'] = annotation_db['comparative'].fillna(False)
    annotation_db = annotation_db[annotation_db['TAREAN_annotation'] != 'rDNA']
    annotation_db = annotation_db.rename({'cluster':'кластер', 'size, %':'размер, %'}, axis=1)
    renamer = {'Putative satellites (low confidence)':'Сателлиты (низковероятные)', 'Putative LTR elements':'Вероятные LTR элементы', 'Putative satellites (high confidence)':'Сателлиты (высоковероятные)'}
    annotation_db['TAREAN_annotation'] = annotation_db['TAREAN_annotation'].apply(lambda x: renamer[x])

    ncbi_df = my_df[my_df['db'].apply(lambda x: 'ncbi' in x)]
    local_df = my_df[my_df['db'].apply(lambda x: ('local' in x) or ('comp' in x))]
    ref_df = my_df[my_df['db'].apply(lambda x: x=='ref')]

    tarean_dict = {}
    for type_gr, group in annotation_db.groupby('TAREAN_annotation'):
        group['размер, %'] = group['размер, %'].apply(lambda x:'{0:.3f}'.format(x))
        group_columns = [x for x in group.columns if not (x in {'Size', 'Size_adjusted', 'TAREAN_annotation', 'comparative', 'Supercluster'} or len(re.findall('^[A-z]+[0-9]+$', x)) > 0)]
        group = group[group_columns]
        ###NEW
#        try:
#            group['gsnap']=group["кластер"].apply(lambda x: glob.glob(f'./igv/{x.replace("(","").replace(")","")}.png')[0])
#            group['pseudo-it']=group["кластер"].apply(lambda x: glob.glob(f'./mapping/igv/{x.replace("(","").replace(")","")}.png')[0])
#        except:
#            pass
        tarean_dict.update({type_gr: group})

    ncbi_dict = {}
    for qseqid, group in ncbi_df.groupby('qseqid'):
        group = group.sort_values(['task', 'fitting', 'sseqid'])[group['sum_coverage']>0.7]
        group['db'] = 'ncbi'
        if len(group) > 0: 
            group = lev(group).head()
        group = group.merge(ncbi_naming, left_on='sseqid', right_on='short').drop('sseqid', axis=1)[['full', 'task', 'sum_coverage']]
        group['sum_coverage'] = group['sum_coverage'].apply(lambda x:'{0:.3f}'.format(x))
        group = group.rename({'full':'повтор из базы NCBI','sum_coverage':'суммарное покрытие'}, axis=1)
        group = group[['повтор из базы NCBI', 'task', 'суммарное покрытие']]
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
        group = group.rename({'sseqid':'повтор из наших сборок','sum_coverage':'суммарное покрытие'}, axis=1)
        group = group[['повтор из наших сборок', 'task', 'суммарное покрытие']]
        found = set(group['повтор из наших сборок'].apply(lambda x: re.findall(' .*_CL', x)[0].replace('_CL','').replace(' ','')).tolist())
        all = set(['KK1', 'KK5(old)', 'KK6(old)', 'KK7(old)'])
        not_found = all-found-set([index, 'KK1KK2KK3KK4'])
        if len(not_found) > 0:
            group.loc[len(group)] = ["*"+",".join(sorted(list(not_found))),'blastn','NOT FOUND']
        local_dict.update({qseqid:group})

    ref_dict = {}
    for qseqid, group in ref_df.groupby('qseqid'):
        group = group.sort_values(['task', 'sseqid'])
        group = group.rename({'sseqid':'олиг из статей или ранее проработанный повтор','sum_coverage':'суммарное покрытие', 'pident':'идентичность', 'length':'длина'}, axis=1)
        group = group[['олиг из статей или ранее проработанный повтор','идентичность', 'длина', 'evalue']]
        ref_dict.update({qseqid:group})

    workbook = xl.Workbook(f'./solo_annotation/{indexi[0]}.xlsx')
    worksheet = workbook.add_worksheet('main_sheet')
    sidecols=2

    get_columns = lambda x: list(x.items())[-1][-1].columns.values.tolist()
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
    type_order = ['Сателлиты (высоковероятные)', 'Сателлиты (низковероятные)', 'Вероятные LTR элементы']
    for key in type_order:
        worksheet.merge_range(row, min(coords.values()), row, max(coords.values())+3, key, merge_format)
        row+=1
        if key in tarean_dict.keys():
            annotation = tarean_dict[key]
        else:
            continue
        for _, an_row in annotation.iterrows():
            cluster = an_row['кластер']
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
                                if 'KK6(old)_' in str(item):
                                    worksheet.write(row+rel_row, col+rel_col, item, workbook.add_format({'bold': True, 'font_color': 'red'}))
                                elif 'KK5(old)_' in str(item):
                                    worksheet.write(row+rel_row, col+rel_col, item, workbook.add_format({'bold': True, 'font_color': 'blue'}))
                                elif 'KK7(old)_' in str(item):
                                    worksheet.write(row+rel_row, col+rel_col, item, workbook.add_format({'bold': True, 'font_color': 'green'}))
                                elif 'KK1_' in str(item):
                                    worksheet.write(row+rel_row, col+rel_col, item.replace('*',''), workbook.add_format({'bold': True, 'font_color': 'purple'}))
                                else:
                                    worksheet.write(row+rel_row, col+rel_col, item)
                            else:
                                if rel_col in [2]: #[2, 5, 6]
                                    path = item
                                    item = ''
                                if row+max(list(map(len, chunks)))-1 != row:
                                    worksheet.merge_range(row, col+rel_col, row+max(list(map(len, chunks)))-1, col+rel_col, item, merge_format)
                                else:
                                    worksheet.write(row, col+rel_col, item, merge_format)
                                if rel_col in [2, ]:#5, 6]:
                                    worksheet.insert_image(row, col+rel_col, path, {'x_scale': 0.15, 'y_scale': 0.15, 'object_position': 1})
            row = row + max(list(map(len, chunks)))

    worksheet.write(1, 0, 'Парных ридов подано на RepEx')
    worksheet.write(2, 0, 'Процент поданных ридов в топе кластеров')
    with open(gl.glob(f'{indexi[0]}/index.html')[0], 'r') as file:
        text = file.read()
        pair_number_context = re.findall("<p class='character'>Number of analyzed reads: [0-9]+</p>", text)[0]
        pair_number_an = int(re.findall('[0-9]+', pair_number_context)[0])
        percent_number_context = re.findall("<p class='character'>Proportion of reads in top clusters : [0-9]+ %</p>", text)[0]
        percent_number = int(re.findall('[0-9]+', percent_number_context)[0])
    worksheet.write(1, 1, pair_number_an)
    worksheet.write(2, 1, percent_number)
    worksheet.merge_range(0,0,0,1, SHDICT[an_row['кластер'].split('_')[0].split('(')[0]], merge_format)

    worksheet.set_column(2,2, 30)
    worksheet.set_column(4,4, 40)
    worksheet.set_column(7,7, 40)
    worksheet.set_column(8,8, 40)
    worksheet.set_column(9,9, 40)
    worksheet.set_column(12,12, 40)
    workbook.close()
