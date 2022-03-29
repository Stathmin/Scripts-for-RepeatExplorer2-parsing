import glob
import os
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('from_path', type=str)
parser.add_argument('to_path', type=str)

args = parser.parse_args()
if not glob.glob(args.to_path):
    os.mkdir(args.to_path)
from_path, to_path = glob.glob(args.from_path)[0], glob.glob(args.to_path)[0]

#files in basic folder
files = [file for x in ['.csv', '.fasta', '.html', '.css'] for file in glob.glob(from_path+'/*'+x)]
basis = [item for item in files if not os.path.isdir(item)]

apply = lambda x, y: list(map(lambda y: x(y), y))
rename = lambda x: x.replace('/tarean', '/graph_layout.png')
annotated_cls = glob.glob(from_path+'/seqclust/clustering/clusters/*/tarean')

#graph layouts for TAREANed repeats
graph_layouts_basic = apply(rename, annotated_cls)
graph_layouts_advanced = [graph
                          if not glob.glob(graph.replace('/graph_layout.png', '/html_files/graph_comparative.png'))
                          else glob.glob(graph.replace('/graph_layout.png', '/html_files/graph_comparative.png'))[0]
                          for graph in graph_layouts_basic 
                          ]


#consensi
get_logo = lambda x: max(glob.glob(x.replace('/tarean', '/tarean/img/logo*')), key=os.path.getctime)
consensi = apply(get_logo, annotated_cls)

ins = basis + graph_layouts_advanced + consensi
rename = lambda x: x.replace('/tarean', '/logo.png')
outs = basis + graph_layouts_basic + apply(rename, annotated_cls)
rename = lambda x: x.replace(from_path, to_path)
outs = apply(rename, outs)


for num, item in enumerate(ins):
    directory = os.path.dirname(outs[num])
    if not os.path.exists(directory):
        os.makedirs(directory)
    shutil.copy2(item, outs[num])