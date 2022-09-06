from ROOT import TMatrixDBase, TMatrixD, TMatrixDSym
import numpy as np

def ArrayToTMatrix(arr):
    shape = np.shape(arr)
    N, M = shape[0], shape[1]
    matr = TMatrixD(N,M)
    for i in range(N):
        for j in range(M):
            matr[i][j] = arr[i][j]
    return matr

def TMatrixToArray(matr):
    N, M = matr.GetNrows(), matr.GetNcols()
    arr = np.zeros([N,M])
    for i in range(N):
        for j in range(M):
            arr[i][j] = matr[i][j]
    return arr

def CalcCov(A_in):
    if isinstance(A_in, TMatrixDBase):
        A = TMatrixToArray(A_in)
    else: A = A_in.copy()
    M = 1./A.shape[1]*np.sum(A,axis=1,keepdims=True)
    B = A-M
    # sample covariance matrix S
    S = 1./(A.shape[1]-1)*np.matmul(B,B.T)
    if isinstance(A_in, TMatrixDBase):
        return ArrayToTMatrix(S)
    else:
        return S

def CovToCorr(cov):
    if isinstance(cov, TMatrixDBase):
        N, M = cov.GetNrows(), cov.GetNcols()
        res = TMatrixDSym(N)
    if isinstance(cov, np.ndarray):
        N, M = cov.shape[0], cov.shape[1]
        res = np.zeros_like(cov)
    assert(N == M)
    for i in range(N):
        for j in range(N):
            if not (cov[i][i]==0 or cov[j][j]==0):
                res[i][j] = cov[i][j]/(np.sqrt(cov[i][i])*np.sqrt(cov[j][j]))
    return res

def TMInvert(matr):
    res = TMatrixD(matr)
    res.Invert()
    return res

def TMTranspose(matr):
    res = TMatrixD(matr)
    res.T()
    return res

def TMSubmatrix(matr,rows=None,cols=None):
    if rows is None: rows=range(matr.GetNrows())
    if cols is None: cols=range(matr.GetNcols())
    res = TMatrixD(len(rows),len(cols))
    for i in range(len(rows)):
        for j in range(len(cols)):
            res[i][j] = matr[rows[i]][cols[j]]
    return res

def TMMultiply(m1,*m2):
    if len(m2)==1:
        res = TMatrixD(m1, TMatrixD.kMult, m2[0])
    elif len(m2)==2:
        res = TMatrixD(m1, TMatrixD.kMult, TMatrixD(m2[0], TMatrixD.kMult, m2[1]))
    return res

def TMTrace(matr):
    assert(matr.GetNrows() == matr.GetNcols())
    N = matr.GetNrows()
    res = 0.
    for i in range(N):
        res += matr[i][i]
    return res
