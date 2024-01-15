import pandas as pd
import requests
import numpy as np
import json 
from sklearn.metrics import jaccard_score

# function to map _id to the corresponding identifier based on _labels
def map_id(row):
    if row['_labels'] == ':Disease':
        return row['diseaseID']
    else:
        return None
    
def jaccard_similarity(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union    

def save_dis_sim(kgfile, sim_file):
    kg_data = pd.read_csv(kgfile, dtype=str)

    # filter rows with ':Disease' in '_labels'
    disease_nodes = kg_data[kg_data['_labels'].isin([':Disease'])].copy()

    # create a new dataframe for neo4j id -> actual id
    id_map = pd.DataFrame({
        'id': disease_nodes['_id'],
        'name': disease_nodes['name'],
        'mapped_id': disease_nodes.apply(map_id, axis=1)
    })

    diseases = id_map["mapped_id"].unique()

    api_host = "https://www.disgenet.org/api"

    api_key = "964e5a12515da18c2b3c4596f71451c68ca0e70a"
    s = requests.Session()

    disease_dict = {}
    for i in range(0, len(diseases)):
        disease_dict[diseases[i]] = []

    if api_key:
        #Add the api key to the requests headers of the requests Session object in order to use the restricted endpoints.
        s.headers.update({"Authorization": "Bearer %s" % api_key})
        print('Disease-gene dict is being generated...')
        query_diseases = ""
        for i in range(0, len(diseases)):
            query_diseases = query_diseases + ',' + diseases[i]
            if i % 100 == 0 or i == len(diseases) - 1:
                gda_response = s.get(api_host+'/gda/disease/' + query_diseases[1:], params={'source':'CURATED'})
                parsed_response = gda_response.json()
                for j in range (0, len(parsed_response)):
                    disease_dict[parsed_response[j]['diseaseid']].append(parsed_response[j]['geneid'])
                query_diseases = ""  
        #print(disease_dict)
        with open('dis_gene_dict.txt','w') as data:  
            data.write(json.dumps(str(disease_dict)))

    if s:
        s.close()

    disease_similarity_matrix = {}

    print('Jaccard similarities are being calculated...')
    for disease1, genes1 in disease_dict.items():
        for disease2, genes2 in disease_dict.items():
            if disease1 != disease2:
                similarity = jaccard_similarity(genes1, genes2)
                disease_similarity_matrix.setdefault(disease1, {})[disease2] = similarity

    sim_graph = pd.read_csv(sim_file, sep="\t")

    print('Similarities are being written...')
    for key, value in disease_similarity_matrix.items():
        for key2, value2 in value.items():
            if value2 > 0:
                sim_graph.loc[len(sim_graph.index)] = [key, key2, 2, value2, len(sim_graph.index)] 

    sim_graph.to_csv(sim_file, sep="\t", index = False, header = False)             
    
