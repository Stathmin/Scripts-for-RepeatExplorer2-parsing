import glob
import re
import os
import subprocess as subprocess
import tempfile
import pandas as pd
import matplotlib as plt
import time
from progress.bar import IncrementalBar
import numpy as np



def last_name(filename):
    #gets file name from path to file, remove format
    return re.sub('\.\w+','', filename.split('/')[-1])

def blast_runner(input_fasta='', db='',
                 blastn_path = 'blastn',
                 task='blastn',
                 num_threads=20,
                 verbose=False,
                 tmpdirname=None):
    #runs blastn task in terminal, stores all files in tmp_folders if given, else at access point
    tmpdirname = str(tmpdirname)+'/' if tmpdirname else ''
    command = f'{blastn_path} -query {input_fasta} -db {db} -task {task} -num_threads {str(num_threads)} -outfmt 6'
    print(command)
    prefix = last_name(input_fasta)+'@'+last_name(db)
    blast_out = open(f'{tmpdirname}{prefix}.out', 'w')
    blast_err = open(f'{tmpdirname}{prefix}.err', 'w')
    try:
        blast_process = subprocess.Popen(command,
                                         shell = True,
                                         stdout = blast_out,
                                         stderr = blast_err)
        blast_process.wait()
    finally:
        if verbose:
            print(command)
        blast_out.close()
        blast_err.close()
    if blast_process.returncode != 0:
        msg = "an error occured during blastn"
        with open(f'{tmpdirname}{prefix}.err', 'r') as blast_err:
            for line in blast_err:
                msg = msg + line
        raise RuntimeError( msg )
    else:
        if verbose:
            print("successful blastn")
        return f'{tmpdirname}{prefix}.out'

def x3_variant(filename, st_dir=os.getcwd()):
    #creates a version of given fasta with trippled seqs at access point or in folder
    with open(filename, 'r') as file:
        txt=file.read()
        txt = '\n'.join(list(map(lambda x: '>'+x.split('\n')[0]+'\n'+(''.join(x.split('\n')[1:]))*3, list(filter(None, txt.split('>'))))))
        x3_path = st_dir + '/' + 'x3_' + last_name(filename) + '.fasta'
        with open(x3_path, 'w') as outfile:
            outfile.write(txt)
        return x3_path

def blast_coordinator(tasks, db_dict, input_fasta, input_fasta_abs, num_threads, verbose, fmt_header):
    #runs blasts at tmp folders with various params, merges their outputs into one pd.DataFrame
    dfs = []
    for db in db_dict.keys():
        for task in tasks:
            with tempfile.TemporaryDirectory() as tmpdirname:
                if '_x3' in db:
                    input_fasta = x3_variant(input_fasta_abs, tmpdirname)
                else:
                    input_fasta = input_fasta_abs
                result_file = blast_runner(input_fasta=input_fasta, db=db_dict[db], task=task, 
                                           verbose=verbose, tmpdirname=tmpdirname, num_threads=num_threads)
                df = pd.read_csv(result_file, names=fmt_header, delimiter='\t')
                df['task'] = task
                df['db'] = db
                dfs.append(df)
    return pd.concat(dfs)


def df_cleaner(df,tasks):
    #drops HSPs with E-value > 0.05, sets column data types, gets lengths FROM NAMES for single and trippled seqs (requires RepeatExpolere naming format for seqs)
    df = pd.DataFrame(df)
    df = df[df['evalue']<=0.01].reset_index().drop('index', axis=1)
    #NEW
    df = df[df['pident']>=90].reset_index().drop('index', axis=1)
    df['task'] = pd.Categorical(df['task'], tasks)
    df['db'] = pd.Categorical(df['db'], ['ref', 'ncbi', 'local', 'comp', 'ncbi_x3', 'local_x3', 'comp_x3'])
    df['qlength'] = df['qseqid'].apply(lambda x: re.findall('\d+', x)[-1]).astype(int)
    df.loc[df['db'].apply(lambda x: '_x3' in x), 'qlength'] = df.loc[df['db'].apply(lambda x: '_x3' in x), 'qlength'].apply(lambda x: x*3)
    df['length'] = df['length'].astype(int)
    return df

def are_connected(chunk1, chunk2):
    #two intervals interection test
    if (min(chunk1) <= min(chunk2) and max(chunk1) >= min(chunk2)):
        return chunk1, chunk2
    elif (min(chunk2) <= min(chunk1) and max(chunk2) >= min(chunk1)):
        return chunk2, chunk1
    else:
        return False

def group_chunks(chunk_list):
    #iterate and merge intersecting intervals while possible
    while True:
        out_list = []
        ignored = [False]*len(chunk_list)
        for num1, chunk1 in enumerate(chunk_list):
            ignored[num1] = True
            for num2, chunk2 in enumerate(chunk_list):
                if not ignored[num2]:
                    intersect = are_connected(chunk1, chunk2)
                    if intersect:
                        out_list.append((min(intersect[0]), max(intersect[1])))
                        ignored[num2] = True
        if out_list == []:
            return chunk_list
        if len(out_list) == len(chunk_list) or len(out_list)==1:
            chunk_list = out_list
            return out_list

        chunk_list = out_list


def sum_lengths(df):
    #calculate cumulative coverage of querry by all HSPs with subject, in bp
    df = df.sort_values('qstart')
    chunks = list(zip(df['qstart'].tolist(), df['qend'].tolist()))
    parts =  group_chunks(chunks)
    summ = sum(list(map(lambda x: int(x[1]) - int(x[0]) + 1, parts)))
    if summ <= 0:
        print('zero!')
    return summ

def return_best_examples(df, evalue_max_olig = 0.001):
    #adds to df new columns with estimated coverage type(), returns the longest interval of df
    total_length = list(set(df['qlength']))[0]
    coverage = max(df['length'])/total_length
    covered_length = sum_lengths(df)
    sum_coverage = covered_length/total_length

    #easy case with chunk length > 0.8
    if coverage > 0.8:
        part = df[df['length'] == max(df['length'])]
        fitting = 'near-full'

    #any chunk > 0.3 and has others, summs up to 0.8
    elif (coverage > 0.3) & (sum_coverage > 0.8):
        part = df[df['length'] == max(df['length'])]
        fitting = 'composite'

    elif max(df['length']) >= 20:
        part = df[df['length'] == max(df['length'])]
        evalue = int(part['evalue'].values.tolist()[0])
        if evalue < evalue_max_olig:
            fitting = 'partial'
        else:
            fitting = 'weak'
    else:
        fitting = 'weak'

    df = df[df['length'] == max(df['length'])]
    df = df.sort_values(['evalue']).head(1)
    out_df = df

    out_df['fitting'] = fitting
    out_df['sum_coverage'] = sum_coverage
    out_df['sum_lengths'] = covered_length
    return out_df


def sort_for_rules(df, subjects):
    #rule 1: near-full > composite > partial > weak
    #rule 2: megablast > dc-megablast > blastn
    #rule 3: x1 > x3
    #rule 4: subjects sort order

    dfs = []
    df['source'] = df['qseqid'].apply(lambda x: x.split('_')[0])
    df['source'] = pd.Categorical(df['source'], subjects)
    df['fitting'] = pd.Categorical(df['fitting'], ['near-full', 'composite', 'partial', 'weak'])
    df['task'] = pd.Categorical(df['task'], ['megablast', 'dc-megablast', 'blastn'])
    df['no_multiplier'] = df['db'].apply(lambda x: '_x3' not in x)
    groups = df.groupby(['qseqid', 'sseqid'])
    for name, table in groups:
        best_item = table.sort_values(['fitting', 'task', 'no_multiplier']).head(1)
        dfs.append(best_item)
    df = pd.concat(dfs).sort_values(['source']).reset_index().drop(['index', 'no_multiplier'], axis=1)
    return df

def levenshtein(seq1, seq2):
    #calculate levenshtein distance between 2 strings
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    matrix = np.zeros ((size_x, size_y))
    for x in range(size_x):
        matrix [x, 0] = x
    for y in range(size_y):
        matrix [0, y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            if seq1[x-1] == seq2[y-1]:
                matrix [x,y] = min(
                    matrix[x-1, y] + 1,
                    matrix[x-1, y-1],
                    matrix[x, y-1] + 1
                )
            else:
                matrix [x,y] = min(
                    matrix[x-1,y] + 1,
                    matrix[x-1,y-1] + 1,
                    matrix[x,y-1] + 1
                )

    return (matrix[size_x - 1, size_y - 1])

def apply_leven(df):
    #drop names with close ncbi accessions to not overflow later tables
    framelist = []
    for cl, frame in df[df['db']=='ncbi'].groupby('qseqid'):
        ncbi_names = sorted(list(set(frame['sseqid'].tolist())))
        list_to_remain=[ncbi_names[0], ]
        for i in ncbi_names[1:]:
            if levenshtein(list_to_remain[-1], i) > 1:
                list_to_remain.append(i)
        list_to_drop = list(set(ncbi_names) - set(list_to_remain))
        frame = frame[~frame['sseqid'].isin(list_to_drop)]
        print(frame)
        framelist.append(frame)
        print(list_to_drop, list_to_remain)
    final_df = pd.concat([df[df['db']!='ncbi'], pd.concat(framelist)])
    return final_df


def filter_by_e_value(df):
    ##select if weak or partial alignment and e-value lower than 0.001
    return df[((df['fitting']!='weak')|(df['fitting']!='partial'))&(df['evalue']<0.001)].reset_index().drop(['index'], axis=1)

def remove_self_blasts(df):
    return df[df['qseqid']!=df['sseqid']].reset_index().drop(['index'], axis=1)

def stairway_view(df):
    encountered_set = set()
    ignore_ids = []
    for index, row in df.iterrows():
        ##not mention later if not from the same object and aligns composite or better at megablast
        if (row['sseqid'].split('_')[0] != row['qseqid'].split('_')[0]) & (((row['fitting']=='near-full') or (row['fitting']=='composite')) & (row['task'] == 'megablast')):
            encountered_set.add(row['sseqid'])
        if row['qseqid'] in encountered_set:
            ignore_ids.append(index)
    return df.iloc[df.index.difference(ignore_ids),:].reset_index().drop(['index'], axis=1)

def parse_fasta(fastapath):
    with open(fastapath,'r') as file:
        txt = file.read()
        txt_list = list(filter(None, txt.split('>')))
        final_dict = {x.split('\n')[0]:'>'+x for x in txt_list}
    return final_dict

def get_from_db(db, item):
    #accesses fasta file, get item by exact name
    parsed_db = parse_fasta(db)
    return parsed_db[item]

def get_from_db_multiple(db, pattern):
    #accesses fasta file, get item by name pattern
    parsed_db = parse_fasta(db)
    patterned_keys = [x for x in list(parsed_db.keys()) if re.match(pattern, x)]
    print(patterned_keys)
    if len(patterned_keys) == 0:
        print(db, pattern)
    return "\n".join((parsed_db[key] for key in patterned_keys))


def protect_brackets(string):
    #needed to find fasta seq.names with regex if brackets in them
    return string.replace('(','\(').replace(')','\)')

def slice_dbs(patterns_list, dbs_list, output_path, db_dict):
    #creates fasta file with given seqs from given dbs
    txt = ''
    for num, pattern in enumerate(patterns_list):
        db = db_dict[dbs_list[num]]
        #print(pattern, db)
        txt += get_from_db_multiple(db, pattern)
    with open(output_path,'w') as file:
        file.write(txt)




def blaster(tasks = ['megablast', 'dc-megablast', 'blastn'], 
            num_threads = 20, 
            verbose = True, 
            subjects = ['KP1', 'KP2', 'KP3', 'KP6', 'KP7', 'KP8', 'KP10'],
            dbs = ['local', 'local', 'local', 'local', 'local', 'local', 'local']
           ):

    db_dict = {
               'ncbi':'./ncbi_repeats_db/ncbi_repeats.fasta',
               'ncbi_x3':'./ncbi_repeats_db/ncbi_repeats_x3.fasta',
               'local':'./local_db_solo/multifasta.fasta',
               'local_x3':'./local_db_solo/multifasta_x3.fasta',
               'comp':'./comparatives_db/COMPBASE.fasta',
               'comp_x3':'./comparatives_db/COMPBASE_x3.fasta',
               'ref':'./important_db/reference.fasta'
               
              }

    patterns = [protect_brackets(x)+'_.*' for x in subjects] # example: KA25_CL1_TR_x1_320nt would be found if KA25 is in subjects

    output_path = 'seqstest.txt'
    slice_dbs(patterns, dbs, output_path, db_dict) #subsets fastas from dbs

    input_fasta = output_path
    input_fasta_abs = os.path.abspath(input_fasta)

    fmt_header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']


    df = blast_coordinator(tasks, db_dict, input_fasta, input_fasta_abs, num_threads, verbose, fmt_header)

    df = df_cleaner(df,tasks)
    print(df)

    dfs = []

    bar = IncrementalBar('Countdown', max = len(df.groupby(['qseqid', 'sseqid', 'task', 'db'])))

    for index, table in df.groupby(['qseqid', 'sseqid', 'task', 'db']):
        bar.next()
        dfs.append(return_best_examples(table, evalue_max_olig=0.001))

    df = pd.concat(dfs).reset_index().drop('index', axis=1)

    bar.finish()

    df = sort_for_rules(df, subjects)
    df = remove_self_blasts(df)
    df = filter_by_e_value(df)

    #df = apply_leven(df)
    if len(subjects) > 1:
        df = stairway_view(df)

    assert (len(df[['qseqid','sseqid']].drop_duplicates()) == len(df[['qseqid','sseqid']]))

    return df
