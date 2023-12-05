import sys
import os
sys.path.append(os.path.abspath('..'))

from DREAMwalk.generate_similarity_net import save_sim_graph
from DREAMwalk.generate_embeddings import save_embedding_files

networkf='../demo/demo_graph.txt'
hierf='../demo/demo_hierarchy.csv'
simf='demo_similarty_graph.txt'
cutoff = 0.4

save_sim_graph(networkf=networkf, hierf=hierf, outputf=simf, cutoff=cutoff)

# node type file should be given for application of heterogeneous Skip-gram
nodetypef='../demo/demo_nodetypes.tsv'
embeddingf='embedding_file.pkl'

save_embedding_files(netf=networkf, sim_netf=simf, outputf=embeddingf,
                    nodetypef=nodetypef, tp_factor=0.3)
