import glob as glob
import pandas as pd
import re as re
import itertools as it
import xlsxwriter as xl
import textwrap as txtw
import numpy as np
import os
import plotly.express as px

from collections import defaultdict

depth_dict=dict((
("Unclassified_repeat (conflicting evidences)",0),
("|--rDNA",1),
("|   |--45S_rDNA",2),
("|   |   |--18S_rDNA",3),
("|   |   |--25S_rDNA",3),
("|   |   '--5.8S_rDNA",3),
("|   '--5S_rDNA",2),
("|--satellite",1),
("'--mobile_element",1),
("|--Class_I",2),
("|   |--SINE",3),
("|   |--LTR",3),
("|   |   |--Ty1_copia",4),
("|   |   |   |--Ale",5),
("|   |   |   |--Alesia",5),
("|   |   |   |--Angela",5),
("|   |   |   |--Bianca",5),
("|   |   |   |--Bryco",5),
("|   |   |   |--Lyco",5),
("|   |   |   |--Gymco-III",5),
("|   |   |   |--Gymco-I",5),
("|   |   |   |--Gymco-II",5),
("|   |   |   |--Ikeros",5),
("|   |   |   |--Ivana",5),
("|   |   |   |--Gymco-IV",5),
("|   |   |   |--Osser",5),
("|   |   |   |--SIRE",5),
("|   |   |   |--TAR",5),
("|   |   |   |--Tork",5),
("|   |   |   '--Ty1-outgroup",5),
("|   |   '--Ty3_gypsy",4),
("|   |       |--non-chromovirus",5),
("|   |       |   |--non-chromo-outgroup",6),
("|   |       |   |--Phygy",6),
("|   |       |   |--Selgy",6),
("|   |       |   '--OTA",6),
("|   |       |       |--Athila",7),
("|   |       |       '--Tat",7),
("|   |       |           |--TatI",8),
("|   |       |           |--TatII",8),
("|   |       |           |--TatIII",8),
("|   |       |           |--Ogre",8),
("|   |       |           '--Retand",8),
("|   |       '--chromovirus",5),
("|   |           |--Chlamyvir",6),
("|   |           |--Tcn1",6),
("|   |           |--chromo-outgroup",6),
("|   |           |--CRM",6),
("|   |           |--Galadriel",6),
("|   |           |--Tekay",6),
("|   |           |--Reina",6),
("|   |           '--chromo-unclass",6),
("|   |--pararetrovirus",3),
("|   |--DIRS",3),
("|   |--Penelope",3),
("|   '--LINE",3),
("'--Class_II",2),
("|--Subclass_1",3),
("|   '--TIR",4),
("|       |--MITE",5),
("|       |--EnSpm_CACTA",5),
("|       |--hAT",5),
("|       |--Kolobok",5),
("|       |--Merlin",5),
("|       |--MuDR_Mutator",5),
("|       |--Novosib",5),
("|       |--P",5),
("|       |--PIF_Harbinger",5),
("|       |--PiggyBac",5),
("|       |--Sola1",5),
("|       |--Sola2",5),
("|       '--Tc1_Mariner",5),
("'--Subclass_2",3),
("'--Helitron",4),
("organelle",0),
("|--plastid",1),
("'--mitochondria",1),
("Unclassified repeat (No evidence)",0),
("contamination",0),
))

#Для Кати

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('index')

args = parser.parse_args()

def get_cyphers():
    with open('cypher.csv', 'r') as file:
        #Шифровки для образцов
        cypherdict = defaultdict(None)
        content = {x.split(',')[0]:x.split(',')[1:] for x in file.read().split('\n')}
        cypherdict.update(content)
        return cypherdict

INDEX=args.index

items = list(map(lambda x: x.split('/')[1], glob.glob(f'./{INDEX}/index.html')))
misc = ['organelle                                       |',
       ' |--plastid                                     |',
       " '--mitochondria                                |",
       'Unclassified repeat (No evidence)               |',
       'contamination                                   |']

outdfs = []
for path, item in zip([glob.glob(f'./{x}/summarized_annotation.html')[0] for x in items], items):
    df = pd.read_fwf(path, skiprows=13, widths=[49, 16, 17, 12, 30], delimiter="\n\t")
    df.columns = ['type', 'percent', 'scl', 'cl', 'reads']
    df = df.drop([74, 75, 79, 80, 82, 83])

    for column in ['percent', 'scl', 'cl', 'reads']:
        df[column] = df[column].apply(lambda x: float(x.replace(' ','').replace('|','').replace('</pre>','')))

    #df['type'] = df['type']#.apply(lambda x: re.sub(r'[ ]+$', '', ''.join(x[0:-1]))) #.replace("'","|"))
    df = df.set_index('type')
    not_clusters = df.loc[misc]['percent'].sum()
    factor = (100 - not_clusters)/100
    df['percent normalized'] = df['percent'].apply(lambda x: round(x*factor,2))
    df.loc[misc,'percent normalized'] = ''

    outdfs.append(df)

df = outdfs[0]
if len(outdfs) > 1:
    for other_df in outdfs[1:]:
        df = pd.merge(left=df, right=other_df, left_on='type', right_on='type')

#df.to_csv('chto.csv')


pd.set_option('display.max_rows', df.shape[0]+1)
for index, row in df.iterrows():
    key = [x for x in list(depth_dict.keys()) if x in index][0]
    df.at[index, 'depth'] = int(depth_dict[key])
df['depth'] = df['depth'].astype(int)

position = ['']*(df['depth'].max()+1)

for index, row in df.iterrows():
    item = row['percent']
    position[row['depth']:] = ['']*(len(position[row['depth']:]))
    position[row['depth']] = index.replace('|','').replace('-','').replace('   ','')
    for ind, it in enumerate(position):
        df.loc[index, f'sum % at level {ind}'] = it
for ind, row in df[[x for x in df.columns if '% at' in x]].iterrows():
    chunk = list(filter(None,row.to_list()))
    type_name = re.sub('^[ ]+','', re.sub('[ ]+$','', chunk[-1])).replace("'",'')
    if len(chunk) > 1:
        parent_name = re.sub('^[ ]+','', re.sub('[ ]+$','', chunk[-2])).replace("'",'')
    else:
        parent_name = ''
    df.at[ind, 'type_name'] = type_name
    df.at[ind, 'parent_name'] = parent_name

for depth_level in reversed([x for x in df.columns if 'sum % at' in x]):
    for item in [x for x in df[depth_level].unique() if x != '']:
        sub_df = df[df[depth_level]==item]
        sub_sum = sub_df['percent normalized'].sum()
        sub_id = sub_df.index[0]
        df.loc[sub_df.index, {depth_level}] = ''
        df.loc[sub_id, {depth_level}] = sub_sum

df['sum percent'] = df[[x for x in df.columns if 'at level' in x]].replace('',0).max(axis=1)
df = df.drop([x for x in df.columns if 'at level' in x], axis=1)
df = df.drop('depth', axis=1)

df.loc[misc,'sum percent'] = ''

import plotly.graph_objects as go

fig =px.sunburst(
    df,
    names='type_name',
    parents='parent_name',
    values='sum percent',
    color='parent_name',
    branchvalues="total",
    hover_data=['percent normalized', 'percent']
)

#fig.update_layout(margin = dict(t=0, l=0, r=0, b=0))

fig.write_image("fig0.png", width=1000, height=1000)

import networkx as nx


Graph_df = nx.from_pandas_edgelist(df, 'type_name', 'parent_name', ['sum percent'], create_using=nx.MultiDiGraph()).reverse()



naming = ' '.join(get_cyphers()[INDEX.replace('(old)','')][0:2])

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter(f'rep_{INDEX}.xlsx', engine='xlsxwriter')

temp_df = df.drop(['parent_name','type_name'],axis=1)
# Convert the dataframe to an XlsxWriter Excel object.
temp_df.to_excel(writer, sheet_name='Sheet1', startrow=1, header=False, startcol=1, index=False)

# Get the xlsxwriter objects from the dataframe writer object.
workbook  = writer.book
worksheet = writer.sheets['Sheet1']

header_format = workbook.add_format({
    'bold': True,
    'text_wrap': False,
    'valign': 'center',
    'align': 'left',
    'border': 1})

for col_num, value in enumerate(temp_df.columns.values):
    worksheet.write(0, col_num + 1, value, header_format)

for row_num, value in enumerate(temp_df.index.values):
    value = re.sub('\s+\|$','',value)
    worksheet.write(row_num + 1, 0, value, header_format)

worksheet.conditional_format('A14:A31', {'type': 'text',
                                     'criteria': 'containing',
                                     'value' : '-', 'format': workbook.add_format({'bg_color': 'green','align': 'left'})})
worksheet.conditional_format('A32:A53', {'type': 'text',
                                     'criteria': 'containing',
                                     'value' : '-', 'format': workbook.add_format({'bg_color': 'yellow','align': 'left'})})
worksheet.conditional_format('A57:G57', {'type': 'text',
                                     'criteria': 'containing',
                                     'value' : '-', 'format': workbook.add_format({'bg_color': 'pink','align': 'left'})})
worksheet.conditional_format('A58:G75', {'type': 'text',
                                     'criteria': 'containing',
                                     'value' : '-', 'format': workbook.add_format({'bg_color': 'blue','align': 'left'})})

worksheet.set_column(0, 0, 40)
worksheet.set_column(5, 6, 15)

worksheet.insert_image('H1','fig0.png')

##SECOND

required_names = ['rDNA',
                 'mobile_element',
                  'Class_I',
                  'Ty1_copia',
                  'Ty3_gypsy',
                  'Class_II',
                  'satellite',
                  'LINE',]
                 #'organelle',
                 #'Unclassified repeat (No evidence)',
                 #'contamination']
#plotted_names= ['LINE', 'rDNA', 'organelle', 'satellite', 'Ty1_copia', 'Ty3_gypsy', 'Class_II', 'Unclassified repeat (No evidence)']
plotted_names= ['LINE', 'rDNA', 'Ty1_copia', 'Ty3_gypsy', 'Class_II']

required_index = [[y for y in df.index if x in y][0] for x in required_names]

renamer = dict(zip(required_index, required_names))

scnd_df = df.loc[required_index, 'sum percent'].rename(renamer, axis=0)
scnd_df.to_excel(writer, sheet_name='Sheet2', startrow=0, header=True, startcol=0, index=True)
worksheet2 = writer.sheets['Sheet2']

fig = px.pie(scnd_df[plotted_names].reset_index(), values='sum percent', names='type', title=f'Типизация ридов {naming} по RepeatExplorer2')
fig.write_image("fig1.png")

worksheet2.insert_image('C4','fig1.png')

optional = []
#for i in ['LINE', 'rDNA', 'organelle', 'contamination']:
for i in ['LINE', 'rDNA']:
    if scnd_df.at[i] != 0:
        if optional:
            break
        else:
            optional.append('Также, найдены')
            break

context_list = ['LINE', 'rDNA', 'organelle', 'contamination']
disc_list = ['элементы типа LINE ({:,.2f}%)', 'rDNA ({:,.2f}%)', 'последовательности органнельной ДНК ({:,.2f}%)', 'неклассифицированные повторяющиеся последовательности ({:,.2f}%)', 'продукты загрязнения адаптером ({:,.2f}%)']
context_dict = dict(zip(context_list, disc_list))

#for num, i in enumerate(['LINE', 'rDNA', 'organelle', 'contamination']):
for num, i in enumerate(['LINE', 'rDNA']):
    if scnd_df.at[i] != 0:
        optional.append(context_dict[i].format(scnd_df.at[i]))


#core=f'''В результате анализа репитома образца {naming} были выявлены следующие типы повторяющихся последовательностей ДНК: {'{:,.2f}'.format(scnd_df.at['mobile_element'])}% мобильных элементов, из которых {'{:,.2f}'.format(scnd_df.at['Class_I'])}% относятся к ретротранспозонам, а {'{:,.2f}'.format(scnd_df.at['Class_II'])}% к ДНК транспозонам. Надсемейство Ty1_copia составило {'{:,.2f}'.format(scnd_df.at['Ty1_copia'])}%, а Ty3_gypsy - {'{:,.2f}'.format(scnd_df.at['Ty3_gypsy'])}%. Доля сателлитов составила {'{:,.2f}'.format(scnd_df.at['satellite'])}%.  Неклассифицированные повторяющиеся последовательности составили {'{:,.2f}'.format(scnd_df.at['Unclassified repeat (No evidence)'])}%'''
core=f'''В результате анализа репитома образца {naming} были выявлены следующие типы повторяющихся последовательностей ДНК: {'{:,.2f}'.format(scnd_df.at['mobile_element'])}% мобильных элементов, из которых {'{:,.2f}'.format(scnd_df.at['Class_I'])}% относятся к ретротранспозонам, а {'{:,.2f}'.format(scnd_df.at['Class_II'])}% к ДНК транспозонам. Надсемейство Ty1_copia составило {'{:,.2f}'.format(scnd_df.at['Ty1_copia'])}%, а Ty3_gypsy - {'{:,.2f}'.format(scnd_df.at['Ty3_gypsy'])}%. Доля сателлитов составила {'{:,.2f}'.format(scnd_df.at['satellite'])}%.'''

txt = core +'\n'+', '.join(optional) + '.'

worksheet2.write('A13',txt)
worksheet2.set_column(0, 2, 20)

##THIRD
out_set=set()
parents = ['Class_I','Ty1_copia','Ty3_gypsy','Class_II']
for parent in parents:
    desc = nx.descendants(Graph_df, parent)
    children_with_value = [x for x in desc if df[df['type_name']==x]['percent'][0]>0]
    out_set.update(children_with_value)

required_names = list(out_set)

required_index = [[y for y in df.index if x in y][0] for x in required_names]

renamer = dict(zip(required_index, required_names))

thrd_df = df.loc[required_index, 'percent'].rename(renamer, axis=0)
thrd_df.to_excel(writer, sheet_name='Sheet3', startrow=0, header=True, startcol=0, index=True)
worksheet3 = writer.sheets['Sheet3']

##TY1
for num, parent in enumerate(parents):
    num = num+2
    desc = nx.descendants(Graph_df, parent)
    children_with_value = [x for x in desc if df[df['type_name']==x]['percent'][0]>0]
    plotted_names = children_with_value
    fig = px.pie(thrd_df[plotted_names].reset_index(), values='percent', names='type', title=f'{naming}, {parent}')
    if not 'Class' in parent:
        fig.write_image(f"fig{num}.png")
        worksheet3.insert_image(f'C{(num-3)*30+4}',f'fig{num}.png')


worksheet3.conditional_format('A1:B7', {'type': 'text',
                                     'criteria': 'containing',
                                     'value' : '-', 'format': workbook.add_format({'bg_color': 'green','align': 'left'})})
worksheet3.conditional_format('A8:B12', {'type': 'text',
                                     'criteria': 'containing',
                                     'value' : '-', 'format': workbook.add_format({'bg_color': 'yellow','align': 'left'})})
worksheet3.conditional_format('A13:B13', {'type': 'text',
                                     'criteria': 'containing',
                                     'value' : '-', 'format': workbook.add_format({'bg_color': 'pink','align': 'left'})})
worksheet3.conditional_format('A14:B16', {'type': 'text',
                                     'criteria': 'containing',
                                     'value' : '-', 'format': workbook.add_format({'bg_color': 'blue','align': 'left'})})

worksheet3.set_column(0, 2, 20)

txt = []

txt.append('Среди мобильных элементов можно выделить два класса, с выделением в первом классе подклассов Ty1_copia и Ty3_gypsy.')
for num,parent in enumerate(['Class_I','Ty1_copia','Ty3_gypsy','Class_II']):
    desc = nx.descendants(Graph_df, parent)
    children_with_value = [x for x in desc if df[df['type_name']==x]['percent'][0]>0]
    if parent =='Class_I':
        children_with_value = set(children_with_value) - set(nx.descendants(Graph_df, 'Ty1_copia')) - set(nx.descendants(Graph_df, 'Ty3_gypsy'))
        children_with_value.update(['Ty1_copia','Ty3_gypsy'])
        children_with_value = list(children_with_value)
    chunks = [f'{x} ({"{:,.2f}".format(df[df["type_name"]==x]["percent"][0])}%)' if (x not in ['Ty1_copia','Ty3_gypsy']) else f'{x} ({"{:,.2f}".format(df[df["type_name"]==x]["sum percent"][0])}%)' for x in children_with_value]
    if df[df["type_name"]==parent]["percent"][0] != 0:
        chunks.append(f'Неидентифицированный {parent} ({"{:,.2f}".format(df[df["type_name"]==parent]["percent"][0])}%)')
    _class = 'Подкласс' if parent in ['Ty1_copia','Ty3_gypsy'] else 'Класс'
    txt.append(f'{_class} {parent} представлен элементами: ' + ', '.join(chunks)+'.')


txt = ' '.join(txt)

worksheet3.write('A20',txt)

writer.save()
