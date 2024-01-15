import numpy as np
import pandas as pd

def generate_files(kg_data):
    ## generate graph and nodetypes files

    # filter rows with ':Protein', ':Drug', or ':Disease' in '_labels'
    nodes_filtered = kg_data[kg_data['_labels'].isin([':Protein', ':Drug', ':Disease'])].copy()
    # create a new dataframe for neo4j id -> actual id
    id_map = pd.DataFrame({
        '_id': nodes_filtered['_id'],
        '_label': nodes_filtered['_labels'],
        'mapped_id': nodes_filtered.apply(map_id, axis=1)
    })

    # create a dictionary mapping _id to mapped_id
    id_map_dict = dict(zip(id_map['_id'], id_map['mapped_id']))
    # create a dictionary mapping mapped_id to _label
    label_map_dict = dict(zip(id_map['mapped_id'], id_map['_label']))

    # filter rows with 'INTERACTS_WITH', 'TARGETS', or 'IS_ASSOCIATED_WITH' in '_type'
    edges_filtered = kg_data[kg_data['_type'].isin(['INTERACTS_WITH', 'TARGETS', 'IS_ASSOCIATED_WITH'])].copy()

    count = 0
    output_graph = pd.DataFrame(columns=['source', 'target', 'edge_type', 'weight', 'edge_id'])
    node_types = {'node': 'type'}

    for index, row in edges_filtered.iterrows():
        if row['_start'] in id_map_dict and row['_end'] in id_map_dict:
            source = id_map_dict.get(row['_start'])
            target = id_map_dict.get(row['_end'])
            edge_type = 1
            if row['_type'] == 'TARGETS':
                edge_type = 2
            elif row['_type'] == 'IS_ASSOCIATED_WITH':
                edge_type = 3
            new_edge_row = {'source': source, 'target': target, 'edge_type': edge_type, 'weight':1, 'edge_id': count}
            output_graph.loc[index] = new_edge_row
            if label_map_dict[source] == ':Protein':
                node_types[source] = 'gene'
            elif label_map_dict[source] == ':Drug':
                node_types[source] = 'drug'
            else:
                node_types[source] = 'disease'
            if label_map_dict[target] == ':Protein':
                node_types[target] = 'gene'
            elif label_map_dict[target] == ':Drug':
                node_types[target] = 'drug'
            else:
                node_types[target] = 'disease'                
            count += 1
    
    # saving the graph as txt file 
    output_graph.to_csv('graph.txt', sep="\t", index = False, header=False)
    print("Graph file is saved!")
    # saving the node types as a tsv file
    with open('nodetypes.tsv', 'w') as f:
        for key in node_types.keys():
            f.write("%s\t%s\n" %(key, node_types[key]))
    print("Node types file is saved!")

    ## generate hierarchy file

    # filter rows with ':Drug' in '_labels'
    drugs_filtered = kg_data[kg_data['_labels'].isin([':Drug'])][['drugID', 'atcClassification']]
    drug_hierarchy = generate_drug_hierarchy(drugs_filtered)
    drug_hierarchy.to_csv('hierarchy.csv', sep=",", index = False)

    # filter rows with ':Drug' in '_labels'
    diseases_filtered = kg_data[kg_data['_labels'].isin([':Disease'])][['diseaseID', 'class']]
    disease_hierarchy = generate_disease_hierarchy(diseases_filtered)

    # writing the merged hierarchy to a csv file
    #hierarchy_frames = [drug_hierarchy, disease_hierarchy]
    #hierarchy_result = pd.concat(hierarchy_frames)
    #hierarchy_result.to_csv('hierarchy.csv', sep=",", index = False)
    print("Hierarchy file is saved!")

    ## generate drug-disease association file

    # filter rows with 'INDICATES' in '_type'
    drug_disease_ints_filtered = kg_data[kg_data['_type'].isin(['INDICATES'])].copy()
    dda_df = pd.DataFrame(columns=['drug','disease','label'])

    for index, row in drug_disease_ints_filtered.iterrows():
        if row['_start'] in id_map_dict and row['_end'] in id_map_dict:
            source = id_map_dict.get(row['_start'])
            target = id_map_dict.get(row['_end'])
            label = 1
            dda_df.loc[len(dda_df) + 1] = {'drug': source, 'disease': target, 'label': label}

    # existing drug_disease pairs
    existing_pairs = set(zip(dda_df['drug'], dda_df['disease']))

    for i in range(1, 11):
    # Create a new DataFrame with equal number of rows from disconnected pairs but with label 0
        new_rows = []
        while len(new_rows) < len(dda_df):
            new_drug = np.random.choice(drugs_filtered['drugID'])
            new_disease = np.random.choice(diseases_filtered['diseaseID'])
            new_pair = (new_drug, new_disease)
            
            if new_pair not in existing_pairs:
                new_rows.append({'drug': new_drug, 'disease': new_disease, 'label': 0})
                # Append the new pair to existing pairs
                existing_pairs.add((new_drug, new_disease))

        # Concatenate the new rows to the original DataFrame
        dda_complete_df = pd.concat([dda_df, pd.DataFrame(new_rows)], ignore_index=True)

        dda_complete_df.to_csv('dda_files/dda' + str(i) + '.tsv', sep="\t", index = False)
        print("Drug-Disease association file is saved!")

# function to map _id to the corresponding identifier based on _labels
def map_id(row):
    if row['_labels'] == ':Protein':
        if pd.isna(row['ncbiID']):
            return row['uniprotID']            
        else:
            return row['ncbiID']
    elif row['_labels'] == ':Drug':
        return row['drugID']
    elif row['_labels'] == ':Disease':
        return row['diseaseID']
    else:
        return None

def generate_drug_hierarchy(drug_df) -> pd.DataFrame:
    drug_hierarchy_df = pd.DataFrame(columns=['child', 'parent'])
    drug_hierarchy_dict = {}
    for index, row in drug_df.iterrows():
        drugID = row['drugID']
        atc_classifications = row['atcClassification'].split(';')
        for atc_classification in atc_classifications:
            atc_classification_list = atc_classification.split(',')
            drug_hierarchy_df.loc[len(drug_hierarchy_df) + 1] = {'child': drugID, 'parent': atc_classification_list[0]}
            for i in range(0, len(atc_classification_list) - 1):
                drug_hierarchy_dict[atc_classification_list[i]] = atc_classification_list[i+1]

    drug_hierarchy_df2 = pd.DataFrame(drug_hierarchy_dict.items(), columns=['child', 'parent'])
    drug_hierarchy_df = pd.concat([drug_hierarchy_df, drug_hierarchy_df2])
    return drug_hierarchy_df

def generate_disease_hierarchy(disease_df) -> pd.DataFrame:
    disease_hierarchy_df = pd.DataFrame(columns=['child', 'parent'])
    disease_hierarchy_dict = {}
    for index, row in disease_df.iterrows():
        diseaseID = row['diseaseID']
        mesh_classifications = row['class'].split(';')
        for mesh_classification in mesh_classifications:
            disease_hierarchy_df.loc[len(disease_hierarchy_df) + 1] = {'child': diseaseID, 'parent': mesh_classification}
            disease_hierarchy_dict[mesh_classification] = mesh_classification[0:1]
            disease_hierarchy_dict[mesh_classification[0:1]] = 'disease'

    disease_hierarchy_df2 = pd.DataFrame(disease_hierarchy_dict.items(), columns=['child', 'parent'])
    disease_hierarchy_df = pd.concat([disease_hierarchy_df, disease_hierarchy_df2])
    return disease_hierarchy_df

def export_files():
    # Read the CSV file into a dataframe
    kg_data = pd.read_csv("covid19-kg.csv", dtype=str)
    print("KG file is loaded!")
    generate_files(kg_data)    