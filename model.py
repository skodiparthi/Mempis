import cobra
import numpy as np



def extract_directions(rxn_list):
    '''
    extract the directions predifined in the reaction model
    
    directions: list
        -1: reverse, +1: forward, 0: reversible
    '''
    directions = np.zeros(len(rxn_list))
    for i, rxn in enumerate(rxn_list):
        if not rxn.reversibility:
            if rxn.lower_bound < 0 and rxn.upper_bound <= 0:
                directions[i] = -1
            else:
                directions[i] = 1
    
    return directions


def reactions_associated_with_genes(sbml_model, gene_list):
    obs_rxn_list = []
    num_knock_outs = 0

    with sbml_model:
        for g in sbml_model.genes:
            if g.id in gene_list:
                pass
            else:
                g.knock_out()
                num_knock_outs += 1
        for i, rxn in enumerate(sbml_model.reactions):
            if rxn.functional:
                # drop out the non-enzymatic reactions (num_genes = 0)
                if len(rxn.genes) > 0:
                    obs_rxn_list.append(rxn.id)
    return obs_rxn_list


def binary_stoich_from_sbml_model(sbml_file, debug=False):
    model = cobra.io.read_sbml_model(sbml_file)
    
    print("num genes:", len(model.genes))
    print("num reactions:", len(model.reactions))
    print("num metabolites:", len(model.metabolites))

    stoich = cobra.util.array.create_stoichiometric_matrix(model)

    ############################################################
    # Create the binary stoichiometric matrix
    ############################################################
    binary_stoich = {}
    all_rxn_list = []
    all_cpd_list = []
    for ridx in range(len(model.reactions)):
        rxn = model.reactions[ridx]
        ## add `rxn_xxxxx_c0` reactions only
        if rxn.id.startswith("rxn") and rxn.id.endswith("c0"):
            if rxn.id not in all_rxn_list:
                all_rxn_list.append(rxn)
            net_stoich = {}
            for m in rxn.metabolites:
                cpd = m.id.split("_")[0]
                if cpd not in all_cpd_list: all_cpd_list.append(cpd)
                if cpd in net_stoich:
                    net_stoich[cpd] += rxn.metabolites[m]
                else:
                    net_stoich[cpd] = rxn.metabolites[m]
            binary_stoich[rxn] = net_stoich

    print("binary_stoich:", len(binary_stoich))
    print("all_rxn_list:", len(all_rxn_list))
    print("all_cpd_list:", len(all_cpd_list))

    # construct a binary stoich
    bin_stoich = np.zeros((len(all_cpd_list), len(all_rxn_list)))
    for rid in binary_stoich:
        net_stoich = binary_stoich[rid]
        cidx = all_rxn_list.index(rid)
        for cpd in net_stoich:
            ridx = all_cpd_list.index(cpd)
            val = net_stoich[cpd]
            if val > 0:
                bin_stoich[ridx, cidx] = 1
            elif val < 0:
                bin_stoich[ridx, cidx] = -1

    # filter out the isolate compounds and null reactions
    num_cpds_in_rxn = np.abs(bin_stoich).sum(axis=0)
    num_rxns_in_cpd = np.abs(bin_stoich).sum(axis=1)

    valid_rxn_idx = np.where(num_cpds_in_rxn!=0)[0]
    valid_cpd_idx = np.where(num_rxns_in_cpd!=0)[0]

    all_rxn_list = np.array(all_rxn_list)
    all_cpd_list = np.array(all_cpd_list)

    rxn_list = all_rxn_list[valid_rxn_idx]
    cpd_list = all_cpd_list[valid_cpd_idx]

    print("rxn_list:", len(rxn_list))
    print("cpd_list:", len(cpd_list))
    
    if debug:
        ## print invalid ones (compounds and reactions only associated with cell transport)
        invalid_rxn_idx = np.where(num_cpds_in_rxn==0)[0]
        invalid_cpd_idx = np.where(num_rxns_in_cpd==0)[0]

        for _rxn in all_rxn_list[invalid_rxn_idx]:
            print(_rxn)

        for _cpd in all_cpd_list[invalid_cpd_idx]:
            print(_cpd)

    stoich = bin_stoich[np.ix_(valid_cpd_idx, valid_rxn_idx)]
    print("Binary stoich:", stoich.shape)
    
    return stoich, list(rxn_list), list(cpd_list)