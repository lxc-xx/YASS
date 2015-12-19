#!/usr/bin/env python
import pprint
import numpy as np
import cPickle as pkl
from scipy.sparse.linalg import arpack
import sys

pp = pprint.PrettyPrinter(depth=6)

# Reference paper: http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-251j-introduction-to-mathematical-programming-fall-2009/readings/MIT6_251JF09_SDP.pdf

# Solve the normal equations for D and y:
# C - theta*X^-1 + theta * X^-1 * D * X^-1 = sum_{i=1}^m y_i A_i
# A_i*D = 0, i = 1,...,m
def solve_normal_equatons(X, A, theta, C):
    # Transfer to MZ=N
    # Solve for Z, where Z:= [D.flatten(), y^T]^T



    m = len(A)
    n = X.shape[0]

    N = np.concatenate((((X - X*C*X/theta).flatten()).T, np.mat(np.zeros((m,1)).astype(float))))

    Q = np.concatenate([ ((-X*a*X/theta).flatten()).T for a in A], axis = 1)

    P = np.concatenate([ a.flatten() for a in A], axis = 0)
    M = np.bmat([[np.mat(np.identity(n*n)), Q],[P, np.mat(np.zeros((m,m)))]])

    #Z = np.linalg.solve(M,N)

    Z, residual, rank, s = np.linalg.lstsq(M,N)
    residual =  np.linalg.norm(N - M*Z)

    D = Z[:n*n, 0].reshape(n,n)
    D = 0.5 * (D+D.T)
    y = Z[n*n:, 0].reshape(m,1)


    #pprint.pprint(A[-2])
    #pprint.pprint(M)
    #pprint.pprint(N)
    #pprint.pprint(Z)
    #pprint.pprint(D)
    #pprint.pprint(y)

    
    return {'D':D,'y':y, 'err': residual}


def check_residual_norm_equations(X,A,theta,C, D, y):
    residual = 0.0

    residual += np.linalg.norm(linear_comb_matrices(y, [X*a*X/theta for a in A]) - X*C*X/theta + X -  D)

    for a in A:
        residual += prod_sum(a,D)

    return residual


def constraint_residual(P):
    b = P['b']
    X = P['X']
    A = P['A']

    b = np.array(b).flatten().tolist()
    res = np.linalg.norm([b_i - prod_sum(X,A_i) for A_i,b_i in zip(A,b)])

    return (res, isPSD(X))

def linear_comb_matrices(y,A):
    yl = y.T.tolist()[0]
    return sum([y_i * a_i for y_i, a_i in zip(yl,A)])
    

def prod_sum(X,Y):
    return np.matrix.sum(np.multiply(X,Y))

def rand_mat(m,n):
    return 2*np.mat(np.random.rand(m,n)) - np.mat(np.ones((m,n)))

def rand_sym_mat(n):
    r = rand_mat(n,n)
    return 0.5*(r+r.T)

def rand_psd_mat(n):
    X = rand_mat(n,n)
    X = X.T*X

    return X

def rand_pd_mat(n):
    X = rand_psd_mat(n)
    X += np.mat(np.eye(n))
    return X

# Generate a random problem
# min C*X
# s.t. A_i*X = b_i, i = 1,...,m; and  X is psd
#
# m: number of inequality constraints
# n: X is nxn matrix
#
# return a dictionary containing A, b, C, and a strictly feasible X

def criterion(C,X):
    return prod_sum(C,X)

def barrier_criterion(C,X,theta):
    return prod_sum(C,X) - theta * np.log(np.linalg.det(X))

def first_deritive_barrier_criterion(C,X,theta):

    return C-theta*X.I

def min_eig(X):
    return np.matrix.min(np.mat(np.linalg.eig(X)[0]))


def isPSD(X):
    try:
        L = np.linalg.cholesky(X)
        return True
    except:
        return False

def trivial_feasible_beta_approx(P, theta = 1.0, beta = 0.25):
    X = P['X'].copy()
    A = P['A']
    C = P['C']
    S = theta * X.I
    L = np.linalg.cholesky(X)
    
    n = X.shape[0]

    stop = np.linalg.norm(np.mat(np.eye(n)) - 1/theta * L.T * S * L)
    #print stop
    

    #solve My=N
    M = np.concatenate([ a.flatten().T for a in A], axis = 1)
    N = (C-theta*X.I).flatten().T

    M , N= M.T*M , M.T*N
    y = np.linalg.solve(M,N)


    P['X'] = X
    P['y'] = y
    P['theta'] = theta
    P['stop'] = stop
    P['beta'] = beta
    P['S'] = S

    return P

def feasible_beta_approx( P , theta = 1.0, beta=0.25):
    X = P['X'].copy()
    A = P['A']
    C = P['C']


    while True:
        L = np.linalg.cholesky(X)
        sol = solve_normal_equatons(X,A,theta,C)

        if sol['err'] > 1e-2: 
            break


        D = sol['D']
        y = np.array(sol['y']).flatten().tolist()

        S = C - sum([y_i*A_i for A_i,y_i in zip(A,y)])

        stop = np.linalg.norm(L.I*D*(L.T).I)

        if stop <= 0.25:
            break
        else:
            alpha = 0.2/stop
            while not isPSD(X+alpha*D):
                alpha *= 0.5
            #print 'barrier_criterion: ' + str(barrier_criterion(C,X,theta))
            P['X'] = X
            #print constraint_residual(P)
            #print 'gap: ' + str(prod_sum(X,S))


            X_prime = X + alpha * D
            X = X_prime
            #print 'stop: ' + str(stop)
            #print 'alpha: ' + str(alpha)

            #t = 1

            ##backtracking line search
            #while True:
            #    #if not isPSD(X-t*D) or barrier_criterion(C,X-t*D,theta) > barrier_criterion(C,X,theta) - 0.5 * t * (np.linalg.norm(first_deritive_barrier_criterion(C,X,theta))**2):
            #    if not isPSD(X-t*alpha*D):
            #        #print isPSD(X-t*alpha*D)
            #        t = 0.25*t
            #    else: 
            #        break

            ##print "stop: " + str(stop)
            ##print "step: " + str(t*alpha)
            ##print "criterion: " + str(criterion(C,X))

            #if t < 1e-4:
            #    break

            #X = X - t * alpha * D

    P['X'] = X
    P['y'] = y
    P['theta'] = theta
    P['stop'] = stop
    P['beta'] = beta
    P['S'] = S
    return P


def sdp(P):
    X = P['X'].copy()
    #S = P['S'].copy()
    A = P['A']
    C = P['C']
    b = P['b']
    theta = P['theta']
    beta = P['beta']
    P['hist'] = []

    n = X.shape[0]
    m = len(A)

    while True:

       #if gap < 1e-4: 
       #    break

       #alpha = 1.0-(np.sqrt(beta) - beta)/(np.sqrt(beta) + np.sqrt(n))
       alpha = 0.8

       theta = alpha * theta
       sol = solve_normal_equatons(X,A,theta,C)

       if sol['err'] > 1e-2: 
           break

       #print 'residual: ' + str(check_residual_norm_equations(X,A,theta,C,sol['D'],sol['y']))

       D = sol['D']
       y = sol['y']
       S = C - linear_comb_matrices(y,A)

       t = 1.0

       while not isPSD(X+t*D):
           t = alpha * t

       #print "theta: " + str(theta)
       primal_criterion = prod_sum(C,X)
       dual_criterion = (np.mat(b)*np.mat(y))[0,0]
       P['hist'].append([primal_criterion, dual_criterion])
       print 'primal criterion: ' + str(prod_sum(C,X))
       print 'dual criterion: ' + str((np.mat(b)*np.mat(y))[0,0])
       print 'duality gap: ' + str(prod_sum(S,X))


       #print 'constraint: ' +  str(constraint_residual(P))
       #print 'step: ' + str(t)

       #gap = prod_sum(X,S) 
       #print 'gap: ' + str(gap)
       #if step < 1e-8:
       #    break

       if theta < 1e-4:
           break 

       X = X + t*D
       #if isPSD(X + D): 
       #    X = X + D
       #    S = theta * X.I
       #else: 
       #    print "oops"
       #    break

    P['X'] = X
    return
    
def test_normal():
    C = np.mat([[2,3],
        [3,2]])
    #A = [np.mat([[1,0],[0,1]]), np.mat([[1,0],[0,2]])]
    A = [np.mat([[1,3],[3,1]]), np.mat([[1,5],[5,1]])]
    X = np.mat([[1,5],[5,1]])
    theta = 0.5

    sol = solve_normal_equatons(X,A,theta,C)

    y = sol['y']
    D = sol['D']

    #print C-theta*X+theta*X.I*D*X.I - linear_comb_matrices(y,A)
    



def random_problem(m,n):
    A = [rand_sym_mat(n) for i in range(m)]

    X = rand_pd_mat(n)

    #X = np.mat(np.diag(np.random.rand(n)))

    b = [prod_sum(X,a) for a in A]
    C = rand_sym_mat(n)

    return {'A':A, 'X':X, 'b':b, 'C':C, 'theta': 1.0, 'beta':0.25}

def parse_problem_file(path):

    mdim = -1
    nblocks = -1
    c = []
    F = []
    ndim = -1
    block_base = []


    with open(path, 'r') as f:
        content = [line.rstrip('\n').strip().split('=')[0] for line in f.readlines()]

        idx = 0
        for line in content:
            if not line.startswith('"'):
                line = line.replace(',',' ')
                if idx == 0:
                    mdim = int(line)
                elif idx ==1:
                    nblock = int(line)
                elif idx ==2:
                    block_struct = [abs(int(x)) for x in line.split()]
                    ndim = sum(block_struct)
                    for i in range(mdim+1):
                        F.append(np.mat(np.zeros((ndim,ndim))))

                    block_base = [sum(block_struct[:i]) for i in range(len(block_struct))]
                elif idx ==3:
                    c = [float(x) for x in line.split()]
                else:
                    vals = line.split()
                    n_mat = int(vals[0])
                    n_block = int(vals[1])
                    rel_i = int(vals[2])
                    rel_j = int(vals[3])
                    val = float(vals[4])

                    abs_i = block_base[n_block-1] + rel_i - 1
                    abs_j = block_base[n_block-1] + rel_j - 1

                    F[n_mat][abs_i,abs_j] = val
                    F[n_mat][abs_j,abs_i] = val

                idx += 1

        C = -F[0]
        b = [-x for x in c]
        A = [-x for x in F[1:]]


        #pp.pprint(C)
        #pp.pprint(b)
        #pp.pprint(A)
        return {'A':A, 'b':b, 'C':C, 'theta': 1.0, 'beta':0.25, 'nblock':nblock, 'block_struct':block_struct, 'block_base':block_base, 'mdim':mdim, 'ndim': ndim}
        


def parse_initial_point(path, P, sparse = True):
    mdim = P['mdim']
    ndim = P['ndim']
    block_struct = P['block_struct']
    block_base = P['block_base']
    nblock = P['nblock']

    with open(path, 'r') as f: 
        content = [line.rstrip('\n').strip().split('=')[0] for line in f.readlines()]

        if sparse:
            idx = 0
            yMat = np.mat(np.zeros((ndim,ndim)))
            xMat = np.mat(np.zeros((ndim,ndim)))
            mats = [xMat, yMat]
            for line in content: 
                line = line.replace(',',' ')
                if not line.startswith('"'):
                    if idx == 0:
                        xVec = [float(x) for x in line.split()]
                    else:
                        vals = line.split()
                        n_mat = int(vals[0])
                        n_block = int(vals[1])
                        rel_i = int(vals[2])
                        rel_j = int(vals[3])
                        val = float(vals[4])

                        abs_i = block_base[n_block-1] + rel_i - 1
                        abs_j = block_base[n_block-1] + rel_j - 1
                        
                        mats[n_mat-1][abs_i,abs_j] = val
                        mats[n_mat-1][abs_j,abs_i] = val
                idx += 1
            P['X'] = yMat
        else:
            xVec_idx = 0
            yMat = np.mat(np.zeros((ndim,ndim)))
            while not content[xVec_idx].startswith('xVec'):
                xVec_idx += 1

            xMat_idx = 0
            while not content[xMat_idx].startswith('xMat'):
                xMat_idx += 1

            yMat_idx = 0
            while not content[yMat_idx].startswith('yMat'):
                yMat_idx += 1

            
            #print content[yMat_idx + 1:]
            #print content[yMat_idx+1:]
            for blk_idx, blk in enumerate(content[yMat_idx + 1:]):
                vals = [float(x) for x in blk.split()]
                blk_size = block_struct[blk_idx]
                base = block_base[blk_idx]

                #print vals
                #print block_struct
                #print blk_idx
                #print blk_size

                for i in range(blk_size):
                    for j in range(blk_size):
                        val = vals[i*blk_size + j]

                        abs_i = base + i
                        abs_j = base + j
                        yMat[abs_i,abs_j] = val


                P['X'] = yMat
    return 0

def main():
    #np.random.seed(41392)
    P = random_problem(50,10)

    print "start C*X: " + str(prod_sum(P['C'],P['X']))

    #P = trivial_feasible_beta_approx(P)
    P = feasible_beta_approx(P)

    sdp(P)

    print "end C*X: " + str(prod_sum(P['C'],P['X']))
    print "If feasible: " + str(constraint_residual(P))

    return 0


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print "Usage: sdp.py problem_file init_file output"
    else:
        P = parse_problem_file(sys.argv[1])
        parse_initial_point(sys.argv[2],P, sparse=False)
        #P = feasible_beta_approx(P)
        sdp(P)
        print "end C*X: " + str(prod_sum(P['C'],P['X']))
        print "If feasible: " + str(constraint_residual(P))

        with open(sys.argv[3], 'w') as f:
            pkl.dump(P,f)



    #main()
    #test_normal()

