from pyomo import environ as pym
import numpy as np
import time
import os


def build_mlp(template_stoich, obs_rxn_idx, obs_cpd_idx, template_intra_met_idx=None, _lambda=3):
    nr = template_stoich.shape[1]
    nx = 2 * nr

    Aref = template_stoich  # original (binary) stoich matrix
    Aaug = np.concatenate((Aref, -Aref), axis=1)  # augmented (binary) stoich matrix 
    print(Aref.shape, Aaug.shape)

    # Constraints

    # 1) "NOT both" forward/backward reactions are activated at the same time
    Aineq = np.zeros((nr, nx))
    bineq = np.zeros(nr)
    for i in range(nr):
        Aineq[i, i] = 1
        Aineq[i, nr+i] = 1
        bineq[i] = 1

    # 2) constrains for measured metabolites
    #    force at least one reaction (around each measured metabolite) to be nonzero
    num_obs_cpds = len(obs_cpd_idx)
    Aineq_ = np.zeros((num_obs_cpds, nx))
    bineq_ = np.zeros(num_obs_cpds)

    for i, cidx in enumerate(obs_cpd_idx):
        nonzero_rxn_idx = np.where(Aaug[cidx, :] != 0)[0]
        Aineq_[i, nonzero_rxn_idx] = -1
        bineq_[i] = -1

    print("Aineq_:", Aineq_.shape, "bineq_:", bineq_.shape)

    Aineq = np.concatenate((Aineq, Aineq_), axis=0)
    bineq = np.concatenate((bineq, bineq_), axis=0)

    print("Aineq:", Aineq.shape, "bineq:", bineq.shape)

    # 3) steady state
    if template_intra_met_idx is None:
        Aeq = Aaug
        beq = np.zeros(Aaug.shape[0])
    else:
        Aeq = Aaug[template_intra_met_idx, :]
        beq = np.zeros(len(template_intra_met_idx))
    
    print("Aeq:", Aeq.shape, "beq:", beq.shape)

    # 4) constrains for measured transcripts (No valid solution)
#     num_obs_rxns = len(obs_rxn_idx)
#     Aeq_ = np.zeros((num_obs_rxns, nx))
#     beq_ = np.ones(num_obs_rxns)
#     for ii in range(num_obs_rxns):
#         Aeq_[ii, obs_rxn_idx[ii]] = 1
#         Aeq_[ii, nr + obs_rxn_idx[ii]] = 1

#     Aeq = np.concatenate((Aeq, Aeq_), axis=0)
#     beq = np.concatenate((beq, beq_), axis=0)

    print("Aeq:", Aeq.shape, "beq:", beq.shape)

    # f = 1 for all x and -lambda for obj_rxn
    f = np.ones(nx)
    for ri in obs_rxn_idx:
        f[ri] -= _lambda
        f[ri+nr] -= _lambda
        
    return Aineq, bineq, Aeq, beq, f



def build_mlp_with_direction(template_stoich, obs_rxn_idx, obs_cpd_idx, rxn_directions, template_intra_met_idx=None, _lambda=2):
    ## Aeq and beq
    Aaug = []
    ridx2A = []
    # store the reaction index and direction for x
    rxns_in_x = []
    
    for ridx in range(len(rxn_directions)):
        direction = rxn_directions[ridx]
        if direction == -1:  # reverse
            ridx2A.append([len(Aaug)])
            Aaug.append(-template_stoich[:,ridx])
            rxns_in_x.append(-ridx)
        elif direction == 1:  # forward
            ridx2A.append([len(Aaug)])
            Aaug.append(template_stoich[:,ridx])
            rxns_in_x.append(ridx)
        else:
            ridx2A.append([len(Aaug), len(Aaug)+1])
            Aaug.append(template_stoich[:,ridx])
            Aaug.append(-template_stoich[:,ridx])
            rxns_in_x.append(ridx)
            rxns_in_x.append(-ridx)

    Aaug = np.array(Aaug).transpose()
    print("Aaug:", Aaug.shape)
    beq = np.zeros(Aaug.shape[0])
    print("beq:", beq.shape)
    
    if template_intra_met_idx is None:
        Aeq = Aaug
    else:
        Aeq = Aaug[template_intra_met_idx, :]
        beq = beq[template_intra_met_idx]

    print("ridx2A:", len(ridx2A))

    Aineq = []
    bineq = []
    for i, a_idx in enumerate(ridx2A):
        if len(a_idx) == 2:
            tmp_Aineq = np.zeros(Aeq.shape[1])
            for j in a_idx:
                tmp_Aineq[j] = 1
            Aineq.append(tmp_Aineq)
            bineq.append(1)
    Aineq = np.array(Aineq)
    bineq = np.array(bineq)
    print("Aineq:", Aineq.shape, "bineq:", bineq.shape)

    num_obs_cpds = len(obs_cpd_idx)
    print("num_obs_cpds:", num_obs_cpds)
    Aineq_ = np.zeros((num_obs_cpds, Aeq.shape[1]))
    bineq_ = np.zeros(num_obs_cpds)
    for i, cidx in enumerate(obs_cpd_idx):
        nonzero_rxn_idx = np.where(Aaug[cidx, :] != 0)[0]
        Aineq_[i, nonzero_rxn_idx] = -1
        bineq_[i] = -1
    print("Aineq_:", Aineq_.shape, "bineq_:", bineq_.shape)

    Aineq = np.concatenate((Aineq, Aineq_), axis=0)
    bineq = np.concatenate((bineq, bineq_), axis=0)
    
    print("Aineq:", Aineq.shape, "bineq:", bineq.shape)

    # reactions associated with over-expressed genes
    # print("reactions associated with over-expressed genes:", len(obs_rxn_list))

    ## Aeq_ and beq_
    # Aeq_ = np.zeros((len(obs_rxn_list), Aeq.shape[1]))
    # beq_ = np.ones(len(obs_rxn_list))
    # valid_rids = [vrxn.id for vrxn in valid_rxn_list]

    # for i, obs_rxn in enumerate(obs_rxn_list):
    #     if obs_rxn in valid_rids:
    #         ridx = valid_rids.index(obs_rxn)
    #         for j in ridx2A[ridx]:
    #             Aeq_[i, j] = 1
    # Aeq = np.concatenate((Aeq, Aeq_), axis=0)
    # beq = np.concatenate((beq, beq_), axis=0)
    print("Aeq:", Aeq.shape, "beq:", beq.shape)

    # f = 1 for all x and -lambda for obj_rxn
    f = np.ones(Aaug.shape[1])
    for ri in obs_rxn_idx:
        for j in ridx2A[ri]:
            f[j] -= _lambda
            
    return Aineq, bineq, Aeq, beq, f, rxns_in_x


def solve(Aineq, bineq, Aeq, beq, f, solver='cplex', email="neos_account_email", solMan="neos"):
    assert f.shape[0] == Aeq.shape[1]
    # provide an email address
    os.environ['NEOS_EMAIL'] = email

    start_time = time.time()

    mip = pym.ConcreteModel()

    mip.I = pym.RangeSet(Aineq.shape[0])
    mip.J = pym.RangeSet(Aeq.shape[0])
    mip.K = pym.RangeSet(Aeq.shape[1])


    mip.x = pym.Var(mip.K, domain=pym.NonNegativeIntegers, bounds=(0,1))

    mip.obj = pym.Objective(expr=pym.quicksum(f[k-1] * mip.x[k] for k in mip.K))


    def ineq_constraint_rule(m, i):
        return pym.quicksum(Aineq[i-1, k-1] * m.x[k] for k in m.K) <= bineq[i-1]

    def eq_constraint_rule(m, j):
        return pym.quicksum(Aeq[j-1, k-1] * m.x[k] for k in m.K) == beq[j-1]

    mip.ConstraintIneq = pym.Constraint(mip.I, rule=ineq_constraint_rule)

    print(">>> ConstraintIneq build time () in {:.2f}s".format(time.time() - start_time))

    mip.ConstraintEq = pym.Constraint(mip.J, rule=eq_constraint_rule)

    print(">>> ConstraintEq build time () in {:.2f}s".format(time.time() - start_time))

    solver_manager = pym.SolverManagerFactory(solMan)
    results = solver_manager.solve(mip, opt=solver)

    print(">>> total time () in {:.2f}s".format(time.time() - start_time))
    
    return [i for i in mip.x if mip.x[i].value!=0]


if __name__ == '__main__':
    FLAGS = parser.parse_args()
    main(FLAGS)