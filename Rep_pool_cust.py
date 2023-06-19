#!/usr/bin/env python3

import argparse as arg
import glob
import os
import subprocess
import sys
#from multiprocessing.pool import ThreadPool
from multiprocessing import Pool
import os


    

def repexjob(INDEX):
    orig= os.getcwd()
    os.chdir(f"/home/dulyanov/RepExTry/SESSIONS_MULT/KP{INDEX}")
    make_command = ["make"]
    stdout = open(f'{INDEX}-STD.out', 'w')
    stderr = open(f'{INDEX}-STD.err', 'w')
    make_proc = subprocess.Popen(make_command, shell = True, stdout = stdout, stderr = stderr)
    make_proc.wait()
    os.chdir(orig)



if __name__=='__main__':
    parser = arg.ArgumentParser()
    parser.add_argument('--adapters', default='/home/ermolaev/soft/bbmap/resources/adapters.fa')
    parser.add_argument('--bbduk_params', default='ktrim=r k=20 mink=10 hdist=2 threads=10 maxns=0 ftl=19 ftr=139 minlen=120')
    parser.add_argument('--num_runs_simult', default='2')
    parser.add_argument('--comparative', action='store_true')
    parser.add_argument('PATHS', help='EXAMPLE: /home/reads/{ind}/{ind}_read_{read_number}.fastq; {ind} and {read_number} should be left unchanged')
    parser.add_argument('INDEXES', nargs='+')
    
    args = parser.parse_args()
    COMPARATIVE, PATHS, ADAPTERS, BBDUK_PARAMS, INDEXES, PROCS = args.comparative, args.PATHS, args.adapters, args.bbduk_params, args.INDEXES, args.num_runs_simult
    

    MAKE_TEMPLATES = glob.glob('/home/dulyanov/RepExTry/**/Makefile_template', recursive=True)
    
    get_path = lambda ind, read_number: ' '.join(sorted(glob.glob(PATHS.format(ind=ind, read_number=read_number))))
    
    if not COMPARATIVE:
        for INDEX in INDEXES:
            SOURCE = locals()
            PATHR1FASTQGZ, PATHR2FASTQGZ = get_path(INDEX, '1'), get_path(INDEX, '2')
            print(PATHS.format(ind=INDEX, read_number='1'))
            REPLACER_DICT = dict((i, SOURCE[i]) for i in ('ADAPTERS', 'BBDUK_PARAMS', 'PATHR1FASTQGZ', 'PATHR2FASTQGZ', 'INDEX'))
            for TEMPLATE in MAKE_TEMPLATES:
                with open(TEMPLATE, 'r') as TFILE:
                    TXT = TFILE.read()
                    for NAME, VALUE in REPLACER_DICT.items():
                        TXT=TXT.replace(f'!{NAME}!', f'{VALUE}')
                RUNNER = TEMPLATE.replace('RepExTry/', f'RepExTry/SESSIONS_MULT/KP{INDEX}/').replace('Makefile_template', 'Makefile')
                os.makedirs(os.path.dirname(RUNNER), exist_ok=True)
                with open(RUNNER, 'w') as RFILE:
                    RFILE.write(TXT)
    else:
        for TEMPLATE in MAKE_TEMPLATES:
            SOURCE = locals()
            with open(TEMPLATE, 'r') as TFILE:
                TXT = TFILE.read()
                if TXT[0:8] == 'ADAPTERS':
                    TXTLIST = TXT.split('\n')
                    BBDUK = TXTLIST[9]
                    TXTLIST.pop(9)
                    for INDEX in INDEXES:
                        SOURCE = locals()
                        PATHR1FASTQGZ, PATHR2FASTQGZ = get_path(INDEX, '1'), get_path(INDEX, '2')
                        REPLACER_DICT = dict((i, SOURCE[i]) for i in ('ADAPTERS', 'BBDUK_PARAMS', 'PATHR1FASTQGZ', 'PATHR2FASTQGZ', 'INDEX'))
                        TXTLIST.insert(9, BBDUK.replace('!PATHR1FASTQGZ!', PATHR1FASTQGZ).replace('!PATHR2FASTQGZ!', PATHR2FASTQGZ).replace('!INDEX!', INDEX))
                    TXT = '\n'.join(TXTLIST)
                for NAME, VALUE in REPLACER_DICT.items():
                    TXT=TXT.replace(f'!{NAME}!', f'{VALUE}')
            RUNNER = TEMPLATE.replace('RepExTry/', f'RepExTry/SESSIONS_MULT/comp{"".join(INDEXES)}/').replace('Makefile_template', 'Makefile')
            os.makedirs(os.path.dirname(RUNNER), exist_ok=True)
    
    pool = Pool(processes=int(PROCS))
    pool.map(repexjob, iter(INDEXES)) 

# with open('Makefile_template', 'r') as file:

