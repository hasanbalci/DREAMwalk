import pandas as pd
import requests
import numpy as np
import json 
from sklearn.metrics import jaccard_score
from itertools import combinations

# function to map _id to the corresponding identifier based on _labels
def map_id(row):
    if row['_labels'] == ':Protein':
        if pd.isna(row['ncbiID']):
            return row['uniprotID']            
        else:
            return row['ncbiID']
    else:
        return None
    
def jaccard_similarity(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union    

def save_protein_sim(kgfile, sim_file):
    kg_data = pd.read_csv(kgfile, dtype=str)

    # filter rows with ':Disease' in '_labels'
    disease_nodes = kg_data[kg_data['_labels'].isin([':Disease'])].copy()

    # filter rows with ':Protein'' in '_labels'
    nodes_filtered = kg_data[kg_data['_labels'].isin([':Protein'])].copy()
    # create a new dataframe for neo4j id -> actual id
    id_map = pd.DataFrame({
        '_id': nodes_filtered['_id'],
        '_label': nodes_filtered['_labels'],
        'mapped_id': nodes_filtered.apply(map_id, axis=1)
    })

    # create a dictionary mapping _id to mapped_id
    id_map_dict = dict(zip(id_map['_id'], id_map['mapped_id']))

    # filter rows with 'IS_INVOLVED_IN' in '_type'
    edges_filtered = kg_data[kg_data['_type'].isin(['IS_INVOLVED_IN'])].copy()

    print('Protein-pathway dict is being generated...')
    protein_dict = {}
    for index, row in edges_filtered.iterrows():
        protein_dict[id_map_dict[row["_start"]]] = []

    for index, row in edges_filtered.iterrows():
        protein_dict[id_map_dict[row["_start"]]].append(row["_end"])

    print('Jaccard similarities are being calculated...')
    protein_similarity_matrix = {}
    protein_combinations = combinations(protein_dict.keys(), 2)

    for protein1, protein2 in protein_combinations:
        similarity = jaccard_similarity(protein_dict[protein1], protein_dict[protein2])
        protein_similarity_matrix.setdefault(protein1, {})[protein2] = similarity 

    print('Similarities are being written...')

    colnames = ['a', 'b', 'c', 'd', 'e'] 
    sim_graph = pd.read_csv(sim_file, names=colnames, sep="\t")
    data = []
    count = len(sim_graph.index)
    for key, value in protein_similarity_matrix.items():
        for key2, value2 in value.items():
            if value2 > 0.4:
                data.append([key, key2, 3, value2, count] )     
                count = count + 1           

    df = pd.DataFrame(data, columns = ['a', 'b', 'c', 'd', 'e']) 
    sim_graph =  pd.concat([sim_graph, df], ignore_index=True, sort=False)

    sim_graph.to_csv(sim_file, sep="\t", index = False, header = False)
    
