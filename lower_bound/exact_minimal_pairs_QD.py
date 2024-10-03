from transition import quasi_differential_transition_matrix
from ortools.sat.python import cp_model
from math import log2, inf
from sage.all import Matrix, vector, GF, span, VectorSpace
from sage.crypto.sbox import SBox
from itertools import product
import numpy as np

S = SBox(0, 2, 4, 12, 8, 14, 24, 21, 16, 19, 28, 5, 17, 20, 11, 23,
         1, 6, 7, 26, 25, 18, 10, 27, 3, 13, 9, 29, 22, 30, 15, 31)
Sinv = S.inverse()
QDTM = quasi_differential_transition_matrix(S, 5, 5)
P = [
    x - 1
    for x in [1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 3, 8, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58, 4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
]
M = Matrix(GF(2), 60, 60)
for i in range(60):
    M[i, i] = 1
    M[i, i-6] = 1
    M[i, i-12] = 1
    M[i, i-19] = 1
    M[i, i-29] = 1
    M[i, i-43] = 1
    M[i, i-51] = 1
Minv = M**-1


def int2vec(x, n=60):
    return vector(GF(2), [(x >> i) & 1 for i in range(n)])


def vec2int(v):
    return sum(int(b) << i for i, b in enumerate(v))


def worst_right_pair_count_subspace(V, d):
    QS = VectorSpace(GF(2), 5)/V
    res = {}
    for x in QS:
        res[tuple(x)] = 0
    for x in range(32):
        x2 = x ^ d
        y = Sinv(x)
        y2 = Sinv(x2)
        if int2vec(y ^ y2, 5) in V:
            res[tuple(QS(int2vec(y, 5)))] += 1
    return min(res.values())

def QD_correlations_subspace(V, d):
    QS = VectorSpace(GF(2), 5)/V
    res = [{} for _ in range(32)]
    for x in QS:
        for s in res:
            s[tuple(x)] = 0
    for x in range(32):
        x2 = x ^ d
        y = Sinv(x)
        y2 = Sinv(x2)
        if int2vec(y ^ y2, 5) in V:
            r = tuple(QS(int2vec(y, 5)))
            for u in range(32):
                res[u][r] += (-1)**((u&x).bit_count())
    lres = tuple(tuple((v, tuple(w for w in x if x[w] == v)) for v in set(x.values())) for x in res)
    return lres

worst_right_pair_count_subspace_lut = {(v, d):worst_right_pair_count_subspace(v, d) for v, d in product(sum((list(span([int2vec(1 << i, 5) for i in range(4)]).subspaces(j)) for j in range(5)), []), range(32))}
QD_correlations_subspace_lut = {(v, d):QD_correlations_subspace(v, d) for v, d in product(sum((list(span([int2vec(1 << i, 5) for i in range(4)]).subspaces(j)) for j in range(5)), []), range(32))}

def optimise_input_space(output_differences, max_dimension, max_subspace): # does not necessarily return all subspaces that have the same maximal minimum number of right pairs over the first sl layer
    dm = max_subspace.dimension()
    model = cp_model.CpModel()
    nb_pairs = 0
    ss_vars = [[model.NewBoolVar("") for _ in range(dm+1)] for _ in output_differences]
    sss = []
    model.Add(sum(sum(i*x for i, x in enumerate(sv)) for sv in ss_vars) <= max_dimension) # dimension restriction
    for i, od in enumerate(output_differences):
        model.Add(sum(ss_vars[i]) == 1)  # only one subspace for each sbox
        sss.append([])
        for dim in range(dm+1):
            mxn = -1
            mxd = []
            for ss in max_subspace.subspaces(dim):
                res = worst_right_pair_count_subspace_lut[(ss, od)]
                if res > mxn:
                    mxn = res
                    mxd = [ss]
                elif res == mxn:
                    mxd.append(ss)
            sss[-1].append(mxd)
            if mxn:
                nb_pairs += log2(mxn)*ss_vars[i][dim]
            else:
                model.Add(ss_vars[i][dim] == 0)
    model.Maximize(nb_pairs)
    solver = cp_model.CpSolver()
    status = solver.solve(model)
    if status == cp_model.OPTIMAL:
        # extract subspaces
        final_subspaces = []
        for i in range(len(output_differences)):
            for dim in range(dm+1):
                if solver.Value(ss_vars[i][dim]):
                    final_subspaces.append(sss[i][dim])
                    break
        return int(round(2**solver.objective_value)), final_subspaces
    else:
        return 0, None

def get_minimal_tweak_pairs_QD(ss, id, od):
    '''
    ss : list of 12 lists of subspaces representing the possible inputs spaces for every sbox of the first sl slayer
    id : list of 12 input differences for the input of the second SL layer
    od: list of 12 output differences for the output of the second SL layer
    returns the minimal number of valid tweak pairs and the input space for every sbox that guarantees this
    will return a negative number of the solver takes more than 5 minutes to find the minimal number of right tweak-pairs
    and will return -inf if the solver can not find a solution within 5 minutes
    '''
    # trails
    input_masks_SL2 = {}
    for masks in product(*[QDTM[32*b, 32*a:32*a+32].nonzero()[0] for a, b in zip(id, od)]):
        c = 1
        for i in range(12):
            c *= QDTM[32*od[i], 32*id[i]+masks[i]]
        input_masks_SL2[masks] = c
    output_masks_SL1 = {}  # map from output masks SL1 to input masks SL2
    for masks in input_masks_SL2:
        v = vector(GF(2), 60)
        for i in range(12):
            v[5*i:5*i+5] = int2vec(masks[i], 5)
        u = v*M
        output_masks_SL1[tuple(vec2int(u[5*i:5*i+5])
                               for i in range(12))] = masks
    b = vector(GF(2), 60)
    for i in range(12):
        b[5*i:5*i+5] = int2vec(id[i], 5)
    a = Minv*b
    od_SL1 = [vec2int(a[5*i:5*i+5]) for i in range(12)]

    # choose input set to minimize the number of trails and simplify the model (this does not guarantee a maximal number of minimal right pairs for the choices of input spaces)
    model = cp_model.CpModel()
    trail_vars = [model.NewBoolVar("") for _ in output_masks_SL1]
    model.Minimize(sum(trail_vars)) # objective, minimize number of trails
    input_space_vars = [[model.NewBoolVar("") for _ in ss[i]] for i in range(12)]
    for i in range(12):
        model.Add(sum(input_space_vars[i]) == 1) # guarantee that only one subspace is selected
    for i, mask in enumerate(output_masks_SL1):
        vars = [model.NewBoolVar("") for _ in range(12)]
        model.AddBoolOr([x.Not() for x in vars] + [trail_vars[i]])
        for j in range(12):
            for k, s in enumerate(ss[j]):
                if len(QD_correlations_subspace_lut[(s, od_SL1[j])][mask[j]]):
                    model.AddImplication(input_space_vars[j][k], vars[j])
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 4
    solver.parameters.max_time_in_seconds = 60
    status = solver.solve(model)
    assert(status == cp_model.OPTIMAL or status == cp_model.FEASIBLE)
    chosen_ss = []
    for i in range(12):
        for j, v in enumerate(input_space_vars[i]):
            if solver.Value(v):
                chosen_ss.append(ss[i][j])
                break

    # find minimal number of right pairs in tweak
    model = cp_model.CpModel()
    key1 = [{tuple(x):model.NewBoolVar("") for x in VectorSpace(GF(2), 5)/chosen_ss[i]} for i in range(12)]
    for i in range(12):
        model.Add(sum(key1[i].values()) == 1)
    key2 = [model.NewBoolVar("") for _ in range(60)]
    acc = 0
    for omasks, imasks2 in output_masks_SL1.items():
        c = input_masks_SL2[imasks2]
        local_pos_values = [QD_correlations_subspace_lut[(chosen_ss[i], od_SL1[i])][omasks[i]] for i in range(12)]
        if any(all(y[0] == 0 for y in x) for x in local_pos_values):
            continue
        ws = []
        w_acc = 0
        for cs in product(*local_pos_values):
            w = model.NewBoolVar("")
            a = 1
            for i, (x, ks) in enumerate(cs):
                # w => ks[i][ks[1]] | ks[i][ks[2]] ...
                model.AddBoolOr([w.Not()] + [key1[i][k] for k in ks])
                a *= x
            w_acc += w
            ws.append((w, a))
        model.Add(w_acc == 1)
        v = model.NewBoolVar("") # u * k2
        model.AddBoolXOr([key2[i] for i in range(60) if (imasks2[i//5]>>(i%5)) & 1] + [v, 1])
        for w, a in ws:
            z = model.NewBoolVar("")
            # z == v * w
            model.AddBoolOr([v.Not(), w.Not(), z])
            model.AddBoolOr([v, z.Not()])
            model.AddBoolOr([w, z.Not()])
            # acc += w*c*(1-2*v)
            acc += a*c*w
            acc -= 2*a*c*z
    model.Minimize(acc)
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = 4
    solver.parameters.max_time_in_seconds = 300
    status = solver.solve(model)
    if status == cp_model.OPTIMAL:
        return solver.objective_value, chosen_ss
    elif status == cp_model.FEASIBLE:
        return -solver.objective_value, chosen_ss
    else:
        return -inf, chosen_ss

def analyze(x, dim=30):
    id_SL2, od_SL2, c = x
    b = vector(GF(2), 60)
    for i in range(12):
        b[5*i:5*i+5] = int2vec(id_SL2[i], 5)
    a = Minv*b
    od_SL1 = tuple(vec2int(a[5*i:5*i+5]) for i in range(12))
    pr0, ss = optimise_input_space(od_SL1, dim, span([int2vec(1 << i, 5) for i in range(4)]))
    pr1 = pr0
    for a, b in zip(id_SL2, od_SL2):
        pr1 *= QDTM[32*b, 32*a]
    if pr1:
        pr2, ss = get_minimal_tweak_pairs_QD(ss, id_SL2, od_SL2)
    else: pr2 = 0
    return id_SL2, od_SL2, c, pr0, pr1, pr2, ss

for d in range(32, 23, -1):
    id_SL2, od_SL2, c, pr0, pr1, pr2, ss = analyze(((0, 0, 10, 0, 9, 0, 0, 0, 0, 0, 0, 0), (0, 0, 4, 0, 16, 0, 0, 0, 0, 0, 0, 0), 6), d)
    print(d, pr0, pr1, pr2)
    print(ss)