import glob
import pandas as pd
import re as re
import itertools as it
import xlsxwriter as xl
import textwrap as txtw
import numpy as np
import os
import re
import argparse
from collections import defaultdict

create_picture = True
try:
    import plotly.express as px
except:
    create_picture = False

#You can select which groups of repeats to exclude from normalization by changing the misc list
misc = ['organelle                                       |',
       ' |--plastid                                     |',
       " '--mitochondria                                |",
       'Unclassified repeat (No evidence)               |',
       'contamination                                   |']

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



parser = argparse.ArgumentParser()
parser.add_argument('index')

args = parser.parse_args()

def get_cyphers():
    with open('cypher.csv', 'r') as file:
        cypherdict = defaultdict(None)
        content = {x.split(',')[0]:x.split(',')[1:] for x in file.read().split('\n')}
        cypherdict.update(content)
        return cypherdict

INDEX=args.index

items = list(map(lambda x: x.split('/')[1], glob.glob(f'./{INDEX}/index.html')))


outdfs = []
for path, item in zip([glob.glob(f'./{x}/summarized_annotation.html')[0] for x in items], items):
    df = pd.read_fwf(path, skiprows=13, widths=[49, 16, 17, 12, 30], delimiter="\n\t")
    df.columns = ['type', 'percent', 'scl', 'cl', 'reads']
    df = df.drop([74, 75, 79, 80, 82, 83])

    for column in ['percent', 'scl', 'cl', 'reads']:
        df[column] = df[column].apply(lambda x: float(x.replace(' ','').replace('|','').replace('</pre>','')))
    df = df.set_index('type')
    nonrepeat_clusters = df.loc[misc]['percent'].sum()
    all_clusters = df['percent'].sum()
    all_non_clusters = 100 - all_clusters
    factor = (100 - all_non_clusters)/100 * (all_clusters - nonrepeat_clusters)/all_clusters
    factor = 1
    df['percent normalized'] = df['percent'].apply(lambda x: round(round(x/factor,6),2))
    df.loc[misc,'percent normalized'] = ''

    outdfs.append(df)

df = outdfs[0]
if len(outdfs) > 1:
    for other_df in outdfs[1:]:
        df = pd.merge(left=df, right=other_df, left_on='type', right_on='type')


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


if create_picture:
    fig =px.sunburst(
        df,
        names='type_name',
        parents='parent_name',
        values='sum percent',
        color='parent_name',
        branchvalues="total",
        hover_data=['percent normalized', 'percent']
    )

    fig.write_html(f"./solo_annotation/{INDEX}_sunburst.html")
    fig.write_image(f"./solo_annotation/{INDEX}.png", width=1000, height=1000)


writer = pd.ExcelWriter(f'./solo_annotation/rep_{INDEX}.xlsx', engine='xlsxwriter')

temp_df = df.drop(['parent_name','type_name'],axis=1)

temp_df.to_excel(writer, sheet_name='Sheet1', startrow=1, header=False, startcol=1, index=False)

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

if create_picture:
    worksheet.insert_image('H1',f"./solo_annotation/{INDEX}.png")


writer.save()