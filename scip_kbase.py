
import numpy as np
import pandas as pd
import glob

############################################################
# DE result files
############################################################
de_files = glob.glob('../SoilSFA_RPKM_Metabolomics_Data/RawCounts_Kbase_*_DE_Result.txt')
oe_genes = []
for de_file in de_files:
    de_df = pd.read_csv(de_file, sep='\t', index_col=0)
    oe_df = de_df[(de_df.padj<=0.05) & (de_df.log2FoldChange>=1)]
    if oe_df.empty:
        print('{}: No Overexpressed genes out of {}'.format(de_file, de_df.shape[0]))
    else:
        print('{}: {} Overexpressed genes out of {}'.format(de_file, oe_df.shape[0], de_df.shape[0]))
        oe_genes += oe_df.index.to_list()
oe_genes = [g.split("_CDS")[0] for g in oe_genes]

############################################################
# metabolites
############################################################
metab_df = pd.read_excel("../SoilSFA_RPKM_Metabolomics_Data/Sums_McClure_MetSuprnt_identifications071321.xlsx", sheet_name='Final')
sample_cols = [c for c in metab_df.columns if "MetSuprnt_Rep" in c]
sample_118h_cols = [c for c in sample_cols if ("MetSuprnt_Rep" in c) & ("_118_hr_" in c)]

cpd_list_obs = metab_df[['cpd'] + sample_118h_cols].dropna().cpd.tolist()
cpd_list_obs = [c + "_c0" for c in cpd_list_obs]

############################################################
# Load the metabolic model from a SBML file
############################################################
import cobra
model = cobra.io.read_sbml_model("../MixedBagComModel.SBML/MixedBagComModel.xml")
print("num genes:", len(model.genes))
print("num reactions:", len(model.reactions))
print("num metabolites:", len(model.metabolites))

stoich = cobra.util.array.create_stoichiometric_matrix(model)

############################################################
## get reactions associated with the over-represented genes
############################################################
reaction_idx = []
num_knock_outs = 0

with model:
    for g in model.genes:
        if g.id in oe_genes:
            pass
        else:
            g.knock_out()
            num_knock_outs += 1
    for i, rxn in enumerate(model.reactions):
        if rxn.functional:
            reaction_idx.append(i)
            

############################################################
## get kbase metabolites
############################################################
observed_midx = []
for midx in range(len(model.metabolites)):
    metabolite = model.metabolites[midx]
    if metabolite.id in cpd_list_obs:
        print(midx, metabolite)
        observed_midx.append(midx)
        
in_metabolites = []
ex_metabolites = []
for midx, metabolite in enumerate(model.metabolites):
    if metabolite.compartment.startswith("e"):
        ex_metabolites.append(midx)
    elif metabolite.compartment.startswith("c"):
        in_metabolites.append(midx)
    else:
        print("compartment is not either `e` or `c`")
print("in_metabolites:", len(in_metabolites))
print("ex_metabolites:", len(ex_metabolites))

######################################################
######################################################
## Aeq and beq
Aeq = []
ridx2A = []
reversible_sets = []  # to get the indices for the reversible pairs
# extract the directions
for ridx in range(len(model.reactions)):
    rxn = model.reactions[ridx]
    if not rxn.reversibility:
        if rxn.lower_bound < 0 and rxn.upper_bound <= 0:
            # directions.append("reverse")
            ridx2A.append([len(Aeq)])
            Aeq.append(-stoich[:,ridx])
        else:
            # directions.append("forward")
            ridx2A.append([len(Aeq)])
            Aeq.append(stoich[:,ridx])
    else:
        reversible_sets.append(len(Aeq))
        ridx2A.append([len(Aeq), len(Aeq)+1])
        Aeq.append(stoich[:,ridx])
        Aeq.append(-stoich[:,ridx])

Aeq = np.array(Aeq).transpose()[np.array(in_metabolites)]
print("Aeq:", Aeq.shape)
print("no. reversible:", len(reversible_sets))
beq = np.zeros(Aeq.shape[0])
print("beq:", beq.shape)

Aineq = np.zeros((len(ridx2A), Aeq.shape[1]))
bineq = np.ones(len(ridx2A))
for i, a_idx in enumerate(ridx2A):
    for j in a_idx:
        Aineq[i, j] = 1
print("Aineq:", Aineq.shape, "bineq:", bineq.shape)

num_obs_cpds = len(observed_midx)
print(len(observed_midx))
Aineq_ = np.zeros((num_obs_cpds, Aeq.shape[1]))
bineq_ = np.zeros(num_obs_cpds)
Aineq_.shape, bineq_.shape
for i, midx in enumerate(observed_midx):
    nonzero_rxn_idx = np.where(Aeq[midx, :] != 0)[0]
    Aineq_[i, nonzero_rxn_idx] = -1
    bineq_[i] = -1
print("Aineq_:", Aineq_.shape, "bineq_:", bineq_.shape)

Aineq = np.concatenate((Aineq, Aineq_), axis=0)
bineq = np.concatenate((bineq, bineq_), axis=0)

print("Aineq:", Aineq.shape, "bineq:", bineq.shape)

######################################################
######################################################

from pyscipopt import Model, quicksum
m = Model("KBase")

num_eq_cons = Aeq.shape[0]
num_ineq_cons = Aineq.shape[0]
num_vars = Aeq.shape[1]

x = {}
for i in range(num_vars):
    name = str(i)
    x[i] = m.addVar(name, vtype='B')

m.setObjective(quicksum(x[k] for k in range(num_vars)))

for i in range(num_eq_cons):
    m.addCons(quicksum(Aeq[i, k] * x[k] for k in range(num_vars)) == beq[i])

for i in range(num_eq_cons):
    m.addCons(quicksum(Aineq[i, k] * x[k] for k in range(num_vars)) <= bineq[i])

m.optimize()
sol = m.getBestSol()

print(sol)

if m.getStatus() != 'optimal':
    print('This is not feasible!')
else:
    print('\nSolution:\n')
    sol = {}
    for i in range(num_vars):
        if m.getVal(x[i]) != 0:
            print(i, m.getVal(x[i]))
        
