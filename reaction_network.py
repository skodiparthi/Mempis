import networkx as nx
import numpy as np


def kegg_reactions_for_graph(kegg_reactions, rxn_idx):
    reactions = {}
    nr = kegg_reactions.shape[0]
    for i in rxn_idx:
        if i > nr:
            ridx = i-1-nr
            direction = -1
        else:
            ridx = i-1
            direction = 1
        rxn = kegg_reactions.iloc[int(ridx)]
        reactions[rxn.rid] = {}
        reactions[rxn.rid]['direction'] = direction
        reactions[rxn.rid]['product'] = rxn['product'].split(",")
        reactions[rxn.rid]['substrate'] = rxn['substrate'].split(",")
    return reactions


def kbase_reactions_for_graph(rxn_list, rxn_idx):
    '''
    rxn_idx: list
        Note: starting from 1 (not from 0)
    '''
    nr = len(rxn_list)
    reactions = {}
    for i in rxn_idx:
        if i > nr:
            ridx = i-1-nr
            direction = -1
        else:
            ridx = i-1
            direction = 1
        rxn = rxn_list[int(ridx)]
        reactions[rxn.id] = {}
        reactions[rxn.id]['direction'] = direction
        reactions[rxn.id]['product'] = [m.id.split("_")[0] for m in rxn.metabolites if np.sign(rxn.metabolites[m])==1 * direction]
        reactions[rxn.id]['substrate'] = [m.id.split("_")[0] for m in rxn.metabolites if np.sign(rxn.metabolites[m])==-1 * direction]
    return reactions


def reactions_for_graph(rxn_list, directions=[]):
    if directions:
        assert len(rxn_list) == len(directions)
    reactions = {}
    for i, rxn in enumerate(rxn_list):
        reactions[rxn.id] = {}
        if directions:
            direction = directions[i]
        else:
            direction = 1
        reactions[rxn.id]['direction'] = direction
        reactions[rxn.id]['product'] = [m.id.split("_")[0] for m in rxn.metabolites if np.sign(rxn.metabolites[m])==1 * direction]
        reactions[rxn.id]['substrate'] = [m.id.split("_")[0] for m in rxn.metabolites if np.sign(rxn.metabolites[m])==-1 * direction]
    return reactions


def create_graph(reactions, rxn_list=[], cpd_list=[]):
    '''
    create the bi-partite graph from the list of reactions
    '''
    
    G = nx.Graph()

    # create nodes (with compounds and reactions)
    all_nodes = []
    for rid in reactions:
        all_nodes.append(rid)
        all_nodes += reactions[rid]['product']
        all_nodes += reactions[rid]['substrate']
    all_nodes = list(set(all_nodes))
    G.add_nodes_from(all_nodes)

    # create edges
    all_edges = []
    for rid in reactions:
        rxn = reactions[rid]
        if rxn['direction'] > 0:
            all_edges += [(rid, cid) for cid in rxn['product']]
            all_edges += [(cid, rid) for cid in rxn['substrate']]
        else:
            all_edges += [(cid, rid) for cid in rxn['product']]
            all_edges += [(rid, cid) for cid in rxn['substrate']]
    G.add_edges_from(all_edges)
    
    ################################################
    count_obs_cpds = []
    count_obs_rxns = []

    color_attributes = {}
    type_attributes = {}

    for node in G.nodes:
        if node.startswith("C") or node.startswith("c"):
            if node in cpd_list:
                count_obs_cpds.append(node)
                color_attributes[node] = {'r': 255, 'g': 0, 'b': 0, 'a': 0.8}
                type_attributes[node] = "Measured Compounds"
            else:
                color_attributes[node] = {'r': 255, 'g': 0, 'b': 0, 'a': 0.3}
                type_attributes[node] = "Added Compounds"
        else:
            if node in rxn_list:
                count_obs_rxns.append(node)
                color_attributes[node] = {'r': 0, 'g': 0, 'b': 255, 'a': 0.8}
                type_attributes[node] = "Measured Reactions"
            else:
                color_attributes[node] = {'r': 0, 'g': 0, 'b': 255, 'a': 0.3}
                type_attributes[node] = "Added Reactions"

    nx.set_node_attributes(G, color_attributes, name="color")

    print("count_obs_cpds:", len(count_obs_cpds))
    print("count_obs_rxns:", len(count_obs_rxns))
    
    return G