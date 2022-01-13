from pymoo.util.dominator import Dominator
import numpy as np
import pandas as pd
import copy

##2. use crowd distance to sort the solution,so each soluton have two new definition 
## caculate crowd distance
## filter_out_duplicates 過濾掉重複項

def fast_non_dominated_sort(F, **kwargs):
    """
    Parameters
    ----------
    F: numpy.ndarray
        objective values for each individual.
    strategy: str
        search strategy, can be "sequential" or "binary".
    Returns
    -------
        fronts: list
            Indices of the individuals in each front.
    References
    ----------
    X. Zhang, Y. Tian, R. Cheng, and Y. Jin,
    An efficient approach to nondominated sorting for evolutionary multiobjective optimization,
    IEEE Transactions on Evolutionary Computation, 2015, 19(2): 201-213.
    copy in pymoo packages by original NSGA-II author 
    """

    
    #M = Dominator.calc_domination_matrix(F)
    M = Dominator.calc_domination_matrix(F.values)

    # calculate the dominance matrix
    n = M.shape[0]

    fronts = []

    if n == 0:
        return fronts

    # final rank that will be returned
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

    # while not all solutions are assigned to a pareto front
    while n_ranked < n:

        next_front = []

        # for each individual in the current front
        for i in current_front:

            # all solutions that are dominated by this individuals
            for j in is_dominating[i]:
                n_dominated[j] -= 1
                if n_dominated[j] == 0:
                    next_front.append(j)
                    ranked[j] = 1.0
                    n_ranked += 1

        fronts.append(next_front)
        current_front = next_front

    return fronts

##2. use crowd distance to sort the solution,so each soluton have two new definition 
## caculate crowd distance
## filter_out_duplicates 過濾掉重複項
def cdist(A, B, **kwargs):
    import scipy
    return scipy.spatial.distance.cdist(A.astype(float), B.astype(float), **kwargs)

def find_duplicates(X, epsilon=1e-16):
    # calculate the distance matrix from each point to another
    D = cdist(X, X)

    # set the diagonal to infinity
    D[np.triu_indices(len(X))] = np.inf

    # set as duplicate if a point is really close to this one
    is_duplicate = np.any(D <= epsilon, axis=1)

    return is_duplicate

def calc_crowding_distance(F, filter_out_duplicates=True):
    n_points, n_obj = F.shape

    if n_points <= 2:
        return np.full(n_points, np.inf)

    else:

        if filter_out_duplicates:
            # filter out solutions which are duplicates - duplicates get a zero finally
            is_unique = np.where(np.logical_not(find_duplicates(F, epsilon=1e-32)))[0]
        else:
            # set every point to be unique without checking it
            is_unique = np.arange(n_points)

        # index the unique points of the array
        #_F = F[is_unique]
        _F = F.values[is_unique]

        # sort each column and get index
        I = np.argsort(_F, axis=0, kind='mergesort')

        # sort the objective space values for the whole matrix
        _F = _F[I, np.arange(n_obj)]

        # calculate the distance from each point to the last and next
        dist = np.row_stack([_F, np.full(n_obj, np.inf)]) - np.row_stack([np.full(n_obj, -np.inf), _F])

        # calculate the norm for each objective - set to NaN if all values are equal
        norm = np.max(_F, axis=0) - np.min(_F, axis=0)
        norm[norm == 0] = np.nan

        # prepare the distance to last and next vectors
        dist_to_last, dist_to_next = dist, np.copy(dist)
        dist_to_last, dist_to_next = dist_to_last[:-1] / norm, dist_to_next[1:] / norm

        # if we divide by zero because all values in one columns are equal replace by none
        dist_to_last[np.isnan(dist_to_last)] = 0.0
        dist_to_next[np.isnan(dist_to_next)] = 0.0

        # sum up the distance to next and last and norm by objectives - also reorder from sorted list
        J = np.argsort(I, axis=0)
        _cd = np.sum(dist_to_last[J, np.arange(n_obj)] + dist_to_next[J, np.arange(n_obj)], axis=1) / n_obj

        # save the final vector which sets the crowding distance for duplicates to zero to be eliminated
        crowding = np.zeros(n_points)
        crowding[is_unique] = _cd

    # crowding[np.isinf(crowding)] = 1e+14
    return crowding

        
#3. Cacaluate Paroto_Optimial Front
def ValuesCal(Total_value,pop_size):
    temp_Paroto_Optimial = fast_non_dominated_sort(Total_value);
    front = np.zeros(pop_size*2,dtype=int);   
    for i in range(len(temp_Paroto_Optimial)):
        for j in temp_Paroto_Optimial[i] :
            front[j] = i;
    fronts = copy.copy(front)
    Total_value["Front_value"] = pd.DataFrame(-front);        
    crowding_of_front = np.zeros(pop_size*2);
    
    for k, fronts in enumerate(temp_Paroto_Optimial):
    # calculate the crowding distance of the front
        crowding_of_front[fronts] = calc_crowding_distance(Total_value.iloc[fronts, :-1])
    crowding_of_front[np.isinf(crowding_of_front)] = -1
##3.1 Compare there front and distenace values  
    Total_value["Crowding_Distance_Value"]= pd.DataFrame(crowding_of_front)    
    return Total_value

#4 Cacaluate Rank 
def RankforNSGAII(Total_value):
    cols = ['Front_value','Crowding_Distance_Value']; #設定要取行名
    tups = Total_value[cols].sort_values(cols,ascending = False).apply(tuple,1) #排名
    f,i = pd.factorize(tups)
    factorized = pd.Series(f+1,tups.index)
    return Total_value.assign(Rank = factorized)
