import warnings

from pymoo.util.dominator import Dominator
from pymoo.core.survival import Survival
from pymoo.docs import parse_doc_string
from pymoo.factory import get_problem, get_reference_directions
import numpy as np
from numpy.linalg import LinAlgError
import pandas as pd
import copy
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting

##use crowd distance to sort the solution,so each soluton have two new definition 

def fast_non_dominated_sort(F, **kwargs):

    
    #M = Dominator.calc_domination_matrix(F)
    M = Dominator.calc_domination_matrix(F.values)

    # calculate the dominance matrix
    n = M.shape[0]

    fronts = []

    if n == 0:
        return fronts

    # the final rank that will be returned
    n_ranked = 0
    ranked = np.zeros(n, dtype=int)

    # for each individual a list of all individuals that are dominated by this one
    is_dominating = [[] for _ in range(n)]

    # storage for the number of solutions dominated this one
    n_dominated = np.zeros(n)

    current_front = []

    for i in range(n):

        for j in range(i + 1, n):
            rel = M[i, j]
            if rel == 1:
                is_dominating[i].append(j)
                n_dominated[j] += 1
            elif rel == -1:
                is_dominating[j].append(i)
                n_dominated[i] += 1

        if n_dominated[i] == 0:
            current_front.append(i)
            ranked[i] = 1.0
            n_ranked += 1

    # append the first front to the current front
    fronts.append(current_front)

    # while not all solutions are assigned to a Pareto front
    while n_ranked < n:

        next_front = []

        # for each individual in the current front
        for i in current_front:

            # all solutions that are dominated by these individuals
            for j in is_dominating[i]:
                n_dominated[j] -= 1
                if n_dominated[j] == 0:
                    next_front.append(j)
                    ranked[j] = 1.0
                    n_ranked += 1

        fronts.append(next_front)
        current_front = next_front

    return fronts


#==============================================================================================================

def compare(a, a_val, b, b_val, method, return_random_if_equal=False):
    if method == 'smaller_is_better':
        if a_val < b_val:
            return a
        elif a_val > b_val:
            return b
        else:
            if return_random_if_equal:
                return np.random.choice([a, b])
            else:
                return None

def calc_perpendicular_distance(N, ref_dirs):
    u = np.tile(ref_dirs, (len(N), 1))
    v = np.tile(N, (len(ref_dirs), 1))

    norm_u = np.linalg.norm(u, axis=1)

    scalar_proj = np.sum(v * u, axis=1) / norm_u
    proj = scalar_proj[:, None] * u / norm_u[:, None]
    val = np.linalg.norm(proj - v, axis=1)
    matrix = np.reshape(val, (len(N), len(ref_dirs)))

    return matrix

def rank_from_fronts(fronts, n):
    rank = np.full(n, 1e16, dtype=int)
    for i, front in enumerate(fronts):
        rank[front] = i

    return rank

# Returns all indices of F that are not dominated by the other objective values
def find_non_dominated(F, _F=None):
    M = Dominator.calc_domination_matrix(F, _F)
    I = np.where(np.all(M >= 0, axis=1))[0]
    return I

def intersect(a, b):
    H = set()
    for entry in b:
        H.add(entry)

    ret = []
    for entry in a:
        if entry in H:
            ret.append(entry)

    return ret

def comp_by_cv_then_random(pop, P, **kwargs):
    S = np.full(P.shape[0], np.nan)

    for i in range(P.shape[0]):
        a, b = P[i, 0], P[i, 1]

        if pop[a].CV > 0.0 or pop[b].CV > 0.0:
            S[i] = compare(a, pop[a].CV, b, pop[b].CV, method='smaller_is_better', return_random_if_equal=True)

        else:
            S[i] = np.random.choice([a, b])

    return S[:, None].astype(int)

def niching(pop, n_remaining, niche_count, niche_of_individuals, dist_to_niche):
    survivors = []

    # boolean array of elements that are considered for each iteration
    mask = np.full(len(pop), True)

    while len(survivors) < n_remaining:

        # number of individuals to select in this iteration
        n_select = n_remaining - len(survivors)

        # all niches where new individuals can be assigned to and the corresponding niche count
        next_niches_list = np.unique(niche_of_individuals[mask])
        next_niche_count = niche_count[next_niches_list]

        # the minimum niche count
        min_niche_count = next_niche_count.min()

        # all niches with the minimum niche count (truncate if randomly if more niches than remaining individuals)
        next_niches = next_niches_list[np.where(next_niche_count == min_niche_count)[0]]
        next_niches = next_niches[np.random.permutation(len(next_niches))[:n_select]]

        for next_niche in next_niches:

            # indices of individuals that are considered and assign to next_niche
            next_ind = np.where(np.logical_and(niche_of_individuals == next_niche, mask))[0]

            # shuffle to break random tie (equal perp. dist) or select randomly
            np.random.shuffle(next_ind)

            if niche_count[next_niche] == 0:
                next_ind = next_ind[np.argmin(dist_to_niche[next_ind])]
            else:
                # already randomized through shuffling
                next_ind = next_ind[0]

            # add the selected individual to the survivors
            mask[next_ind] = False
            survivors.append(int(next_ind))

            # increase the corresponding niche count
            niche_count[next_niche] += 1

    return survivors


def associate_to_niches(F, niches, ideal_point, nadir_point, utopian_epsilon=0.0):
    utopian_point = ideal_point - utopian_epsilon

    denom = nadir_point - utopian_point
    denom[denom == 0] = 1e-12

    # normalize by ideal point and intercepts
    N = (F - utopian_point) / denom
    dist_matrix = calc_perpendicular_distance(N, niches)

    niche_of_individuals = np.argmin(dist_matrix, axis=1)
    dist_to_niche = dist_matrix[np.arange(F.shape[0]), niche_of_individuals]

    return niche_of_individuals, dist_to_niche, dist_matrix


def calc_niche_count(n_niches, niche_of_individuals):
    niche_count = np.zeros(n_niches, dtype=int)
    index, count = np.unique(niche_of_individuals, return_counts=True)
    niche_count[index] = count
    return niche_count


class HyperplaneNormalization:

    def __init__(self, n_dim) -> None:
        super().__init__()
        self.ideal_point = np.full(n_dim, np.inf)
        self.worst_point = np.full(n_dim, -np.inf)
        self.nadir_point = None
        self.extreme_points = None

    def update(self, F, nds=None):

        # find or usually update the new ideal point - from feasible solutions
        self.ideal_point = np.min(np.vstack((self.ideal_point, F)), axis=0)
        self.worst_point = np.max(np.vstack((self.worst_point, F)), axis=0)

        # this decides whether only non-dominated points or all points are used to determine the extreme points
        if nds is None:
            nds = np.arange(len(F))

        # find the extreme points for normalization
        self.extreme_points = get_extreme_points_c(F[nds, :], self.ideal_point,
                                                   extreme_points=self.extreme_points)

        # find the intercepts for normalization and do backup if gaussian elimination fails
        worst_of_population = np.max(F, axis=0)
        worst_of_front = np.max(F[nds, :], axis=0)

        self.nadir_point = get_nadir_point(self.extreme_points, self.ideal_point, self.worst_point,
                                           worst_of_population, worst_of_front)


def get_extreme_points_c(F, ideal_point, extreme_points=None):
    # calculate the asf which is used for the extreme point decomposition
    weights = np.eye(F.shape[1])
    weights[weights == 0] = 1e6

    # add the old extreme points to never loose them for normalization
    _F = F
    if extreme_points is not None:
        _F = np.concatenate([extreme_points, _F], axis=0)

    # use __F because we substitute small values to be 0
    __F = _F - ideal_point
    __F[__F < 1e-3] = 0

    # update the extreme points for the normalization having the highest asf value each
    F_asf = np.max(__F * weights[:, None, :], axis=2)

    I = np.argmin(F_asf, axis=1)
    extreme_points = _F[I, :]

    return extreme_points


def get_nadir_point(extreme_points, ideal_point, worst_point, worst_of_front, worst_of_population):
    try:

        # find the intercepts using gaussian elimination
        M = extreme_points - ideal_point
        b = np.ones(extreme_points.shape[1])
        plane = np.linalg.solve(M, b)

        warnings.simplefilter("ignore")
        intercepts = 1 / plane

        nadir_point = ideal_point + intercepts

        # check if the hyperplane makes sense
        if not np.allclose(np.dot(M, plane), b) or np.any(intercepts <= 1e-6):
            raise LinAlgError()

        # if the nadir point should be larger than any value discovered so far set it to that value
        b = nadir_point > worst_point
        nadir_point[b] = worst_point[b]

    except LinAlgError:

        # fall back to worst of front otherwise
        nadir_point = worst_of_front

    # if the range is too small set it to worst of population
    b = nadir_point - ideal_point <= 1e-6
    nadir_point[b] = worst_of_population[b]

    return nadir_point


#============================================================================================================
#Cacaluate Paroto_Optimial Front
def ValuesCal(Total_value,pop_size):
    temp_Paroto_Optimial = fast_non_dominated_sort(Total_value);
    front = np.zeros(pop_size,dtype=int);
    for i in range(len(temp_Paroto_Optimial)):
        for j in temp_Paroto_Optimial[i] :
            front[j] = i;
    fronts = copy.copy(front)
    
    for i in range(len(fast_non_dominated_sort(Total_value))):
        for j in range(len(fast_non_dominated_sort(Total_value)[i])):
            fronts[fast_non_dominated_sort(Total_value)[i][j]] = i
    Total_value["Front_value"] = pd.DataFrame(fronts);
    crowding_of_front = np.zeros(pop_size);
    #==========================================================================================================
    ref_dirs = get_reference_directions("das-dennis", 2, n_partitions=99)
    global niches
    global pop
    global n_survive
    # attributes to be set after the survival
    F = Total_value["Front_value"]

    # calculate the fronts of the population
    test_fast_non = fast_non_dominated_sort(Total_value)
    non_dominated, last_front = test_fast_non[0], test_fast_non[-1]

    # update the hyperplane based boundary estimation
    hyp_norm = HyperplaneNormalization(ref_dirs.shape[1])
    hyp_norm.update(np.array(Total_value.iloc[:,:-1]), nds=np.array(non_dominated))
    ideal, nadir = hyp_norm.ideal_point, hyp_norm.nadir_point

    # consider only the population until we come to the splitting front
    I = np.concatenate(test_fast_non)
    Total_value, F = Total_value.sort_values(by=['Front_value']), F[I]

    # update the front indices for the current population
    counter = 0
    for i in range(len(test_fast_non)):
        for j in range(len(test_fast_non[i])):
            test_fast_non[i][j] = counter
            counter += 1
    last_front = test_fast_non[-1]

    # associate individuals to niches
    niche_of_individuals, dist_to_niche, dist_matrix = \
        associate_to_niches(Total_value.iloc[:,:-1], ref_dirs, ideal, nadir)

    # attributes of a population
    Total_value.insert(3,'niche', niche_of_individuals)
    Total_value.insert(4,'dist_to_niche', dist_to_niche)

    # set the optimum, first front and closest to all reference directions
    closest = np.unique(dist_matrix[:, np.unique(niche_of_individuals)].argmin(axis=0))
    set_optimum = Total_value.iloc[intersect(test_fast_non[0], closest),:]
    
    Total_value['Rank'] = pd.DataFrame([0]*pop_size)
    for k in range(0,pop_size):
        Total_value.iloc[k,5] = k+1
    
    return Total_value