from urllib.request import urlopen
import pandas as pd
import json
import scipy as sp
import scipy.stats as stats

def get_ko_to_ec():
    '''
    get mapping table from KO to EC
    '''
    return get_link('ec', 'ko')

def get_ko_to_rn():
    '''
    get mapping table from KO to reactions
    '''
    return get_link('rn', 'ko')
    
    
def get_cpd_to_rn():
    '''
    get mapping table from compounds to reactions
    '''
    return get_link('rn', 'cpd')
    
    
def get_cpd_to_ec():
    '''
    get mapping table from compounds to EC
    '''
    return get_link('ec', 'cpd')
    
    
def get_link(target, source):
    '''
    get a mapping table from source to target dbs
    '''
    link = "http://rest.kegg.jp/link/{}/{}".format(target, source)
    source_to_target = []
    
    with urlopen(link) as f:
        for l in f:
            tk = l.decode("utf-8").strip().split("\t")
            source_to_target.append({source: tk[0], target: tk[1]})
    return pd.DataFrame(source_to_target)

def get_pathways():
    '''
    get kegg pathways
    '''
    link = "http://rest.kegg.jp/list/pathway"
    kegg_pathways = []
    with urlopen(link) as f:
        for l in f:
            tk = l.decode("utf-8").strip().split("\t")
            kegg_pathways.append({"pid": tk[0], "desc": tk[1]})
    return pd.DataFrame(kegg_pathways)

def get_brite_json(brite_id):
    '''
    get a brite hierarchy file
    '''
    link = "http://rest.kegg.jp/get/br:{}/json".format(brite_id)
    kegg_pathways = []
    try:
        with urlopen(link) as f:
            return json.load(f)
    except Exception as e:
        return None
    
def get_pathway_with_hierarchy():
    '''
    get a kegg pathway hierarchy
    '''
    data = get_brite_json('br08901')
    if data:
        pathway_hierarchy = []
        for i in data['children']:
            # print("-+"*50)
            # print(i['name'])
            for j in i['children']:
                # print("    ", j['name'])
                for k in j['children']:
                    # print("         ", k)
                    pid, desc = k['name'].split("  ")
                    pathway_hierarchy.append({"Class": i['name'], "Subclass": j['name'], "pid": "path:map"+pid, "desc": desc})
        return pd.DataFrame(pathway_hierarchy)
    return None

def get_ko_numbers():
    '''
    get KEGG Orthology (KO) numbers
    '''
    link = "http://rest.kegg.jp/list/ko"
    kegg_kos = []
    with urlopen(link) as f:
        for l in f:
            tk = l.decode("utf-8").strip().split("\t")
            kegg_kos.append({"ko": tk[0], "desc": tk[1]})
    return pd.DataFrame(kegg_kos)

def get_ko2pathway():
    '''
    get a mapping table between KO numbers and KEGG pathways
    '''
    return get_link('pathway', 'ko')

def get_rxn2pathway():
    '''
    get a mapping table between KEGG reactions and pathways
    '''
    return get_link('pathway', 'rn')

def kegg_fisher_test(ko_list, target_pid, total_num_kos, ko2pathway, alternative='greater'):
    total_num_obs_kos = len(ko_list)
    
    a = ko2pathway[(ko2pathway.pid==target_pid) & (ko2pathway.ko.isin(["ko:"+ko for ko in ko_list]))].shape[0]
    b = total_num_obs_kos - a
    c = ko2pathway[(ko2pathway.pid==target_pid)].shape[0] - a
    d = total_num_kos - (a + b + c)
    
    return stats.fisher_exact([[a, b], [c, d]], alternative=alternative), (a,b,c,d)