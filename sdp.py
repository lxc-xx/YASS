#!/usr/bin/env python
import numpy as np
from scipy.sparse.linalg import arpack

# Reference paper: http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-251j-introduction-to-mathematical-programming-fall-2009/readings/MIT6_251JF09_SDP.pdf

# Solve the normal equations for D and y:
# C - theta*X^-1 + theta * X^-1 * D * X^-1 = sum_{i=1}^m y_i A_i
# A_i*D = 0, i = 1,...,m
def solve_normal_equatons(X, A, theta, C):
    # Transfer to MZ=N
    # Solve for Z, where Z:= [D.flatten(), y^T]^T

    m = len(A)
    n = X.shape[0]

    N = np.concatenate((((X - X*C*X).flatten()).T, np.mat(np.zeros((m,1)).astype(float))))
    Q = np.concatenate([ ((-X*a*X).flatten()).T for a in A], axis = 1)
    P = np.concatenate([ a.flatten() for a in A], axis = 0)
    M = np.bmat([[np.mat(np.identity(n*n)), Q],[P, np.mat(np.zeros((m,m)))]])

    #print M.shape
    #print N.shape
    Z = np.linalg.solve(M,N)

    D = Z[:n*n, 0].reshape(n,n)
    y = Z[n*n:, 0].reshape(m,1)

    return {'D':D,'y':y}


def constraint_residual(P):
    b = P['b']
    X = P['X']
    A = P['A']

    b = np.array(b).flatten().tolist()
    res = np.linalg.norm([b_i - prod_sum(X,A_i) for A_i,b_i in zip(A,b)])

    return (res, isPSD(X))

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
def random_problem(m,n):
    A = [rand_sym_mat(n) for i in range(m)]

    #X = rand_pd_mat(n)

    X = 100*np.mat(np.diag(np.random.rand(n)))

    b = [prod_sum(X,a) for a in A]
    C = rand_sym_mat(n)

    return {'A':A, 'X':X, 'b':b, 'C':C}

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
        D = sol['D']
        y = np.array(sol['y']).flatten().tolist()
        S = C - sum([y_i*A_i for A_i,y_i in zip(A,y)])
        stop = np.linalg.norm(L.I*D*(L.I).T)

        if stop <= 0.25:
            break
        else:
            alpha = 0.2/stop
            t = 1

            #backtracking line search
            while True:
                #if not isPSD(X-t*D) or barrier_criterion(C,X-t*D,theta) > barrier_criterion(C,X,theta) - 0.5 * t * (np.linalg.norm(first_deritive_barrier_criterion(C,X,theta))**2):
                if not isPSD(X-t*alpha*D):
                    print isPSD(X-t*alpha*D)
                    t = 0.25*t
                else: 
                    break

            print "stop: " + str(stop)
            print "step: " + str(t*alpha)
            print "criterion: " + str(criterion(C,X))

            if t < 1e-4:
                break

            X = X - t * alpha * D

    P['X'] = X
    P['y'] = y
    P['theta'] = theta
    P['stop'] = stop
    P['beta'] = beta
    P['S'] = S
    return P


def sdp(P):
    X = P['X'].copy()
    S = P['S'].copy()
    A = P['A']
    C = P['C']
    b = P['b']
    theta = P['theta']
    beta = P['beta']

    n = X.shape[0]
    m = len(A)

    while True:
       gap = prod_sum(X,S) 
       print 'gap: ' + str(gap)
       print 'criterion: ' + str(prod_sum(C,X))

       if gap < 1e-4: 
           break

       alpha = 1.0-(np.sqrt(beta) - beta)/(np.sqrt(beta) + np.sqrt(n))

       theta = alpha * theta
       print "theta: " + str(theta)
       sol = solve_normal_equatons(X,A,theta,C)

       D = sol['D']
       y = np.array(sol['y']).flatten().tolist()
       #S = C - sum([y_i*A_i for A_i,y_i in zip(A,y)])

       #if np.linalg.norm(D) < 1e-4:
       #    break

       if True: #isPSD(X + D): 
           X = X + D
           X = 0.5*(X+X.T)
           S = theta * X.I
       else: 
           print "oops"
           break

    P['X'] = X
    return
    
def main():
    #np.random.seed(100402)
    P = random_problem(5,5)

    print "start C*X: " + str(prod_sum(P['C'],P['X']))
    P = trivial_feasible_beta_approx(P)
    sdp(P)
    print "end C*X: " + str(prod_sum(P['C'],P['X']))
    print constraint_residual(P)

    return 0

if __name__ == "__main__":
    main()

