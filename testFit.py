import yaml
import numpy as np
import ROOT
import sys
from math import sqrt


def ReadIndependent(entry, col=0):
    values = entry['independent_variables'][col]['values']
    if 'value' in values[0]:
        return [X['value'] for X in values]
    else:
        return [(X['low'], X['high']) for X in values]


def ReadDependent(entry, col=0, error=list(), sym_errors=True):
    vals = entry['dependent_variables'][col]['values']
    if len(error) > 0:
        res = list()
        for v in vals:
            sum_hi_sq = 0.
            sum_lo_sq = 0.
            for ecol in error:
                err = v['errors'][ecol]
                if 'symerror' in err:
                    sum_hi_sq += pow(err['symerror'], 2)
                    sum_lo_sq += pow(err['symerror'], 2)
                elif 'asymerror' in err:
                    sum_hi_sq += pow(err['asymerror']['plus'], 2)
                    sum_lo_sq += pow(err['asymerror']['minus'], 2)
            sum_hi = sqrt(sum_hi_sq)
            sum_lo = sqrt(sum_lo_sq)
            if sym_errors:
                res.append((sum_lo + sum_hi) / 2.)
            else:
                res.append((-1. * sum_lo, +1. * sum_hi))
        return np.array(res)
    else:
        return np.array([X['value'] for X in entry['dependent_variables'][col]['values']])


def ReadHgg():
    with open("hepdata_inputs/HEPData-ins1851456-v2-STXS_stage_1.2_minimal_merging_scheme.yaml", "r") as f:
        vals = yaml.safe_load(f)
    with open("hepdata_inputs/HEPData-ins1851456-v2-Correlations__STXS_stage_1.2_minimal_merging_scheme.yaml", "r") as f:
        corr = yaml.safe_load(f)

    labels = ReadIndependent(vals, col=0)
    N = len(labels)
    sm = ReadDependent(vals, col=0)
    bf = ReadDependent(vals, col=1)
    bf_unc = ReadDependent(vals, col=1, error=[0], sym_errors=True)
    assert(N == len(sm) == len(bf) == len(bf_unc))

    bf = bf / sm
    bf_unc = bf_unc / sm

    corr_vals = ReadDependent(corr, col=0)
    assert((N * N) == len(corr_vals))

    cor = ROOT.TMatrixDSym(N)
    cov = ROOT.TMatrixDSym(N)

    for i in range(N):
        for j in range(N):
            cor[i][j] = corr_vals[i * N + j]
            cov[i][j] = corr_vals[i * N + j] * bf_unc[i] * bf_unc[j]
    return {
        'N': N,
        'labels': labels,
        'sm': sm,
        'bf': bf,
        'bf_unc': bf_unc,
        'cor': cor,
        'cov': cov
    }


def ReadWg():
    def Reshape(arr):
        return arr.reshape((3, int(len(arr) / 3))).T.reshape((-1))

    with open("hepdata_inputs/HEPData-ins1978840-v1-p_{mathrm{T}}^{gamma}_times_|phi_{f}|_cross_section.yaml", "r") as f:
        vals = yaml.safe_load(f)
    with open("hepdata_inputs/HEPData-ins1978840-v1-p_{mathrm{T}}^{gamma}_times_|phi_{f}|_correlations.yaml", "r") as f:
        corr = yaml.safe_load(f)

    # labels_j = ReadIndependent(vals, col=0)
    labels_i = ReadIndependent(vals, col=1)
    labels = list()
    N = len(labels_i)
    for i in range(int(N / 3)):
        for j in range(3):
            labels.append('WG_pt_%i_phi_%i' % (i, j))
    # reorder to give pt first, then phi
    sm = Reshape(ReadDependent(vals, col=1))
    bf = Reshape(ReadDependent(vals, col=0))
    bf_unc = Reshape(ReadDependent(vals, col=0, error=[0, 1], sym_errors=True))
    assert(N == len(sm) == len(bf) == len(bf_unc))

    bf = bf / sm
    bf_unc = bf_unc / sm

    corr_vals = ReadDependent(corr, col=0)
    # assert((N * N) == len(corr_vals))

    cor = ROOT.TMatrixDSym(N)
    cov = ROOT.TMatrixDSym(N)

    idx = 0
    for i in range(N):
        for j in range(N - i):
            print(i, N - j - 1, corr_vals[idx])
            cor[i][N - j - 1] = corr_vals[idx]
            cor[N - j - 1][i] = cor[i][N - j - 1] 
            cov[i][N - j - 1] = corr_vals[idx] * bf_unc[i] * bf_unc[N - j - 1]
            cov[N - j - 1][i] = cov[i][N - j - 1]
            idx += 1

    return {
        'N': N,
        'labels': labels,
        'sm': sm,
        'bf': bf,
        'bf_unc': bf_unc,
        'cor': cor,
        'cov': cov
    }


def ReadSingleT():
    with open("hepdata_inputs/HEPData-ins1744604-v1-Table_1.yaml", "r") as f:
        vals = yaml.safe_load(f)
    with open("hepdata_inputs/HEPData-ins1744604-v1-Table_2.yaml", "r") as f:
        corr = yaml.safe_load(f)

    labels = ReadIndependent(vals, col=0)
    labels = ['pt_t_%i_%i' % X for X in labels]
    N = len(labels)
    # sm = ReadDependent(vals, col=0)
    bf = ReadDependent(vals, col=0)
    # TODO: analysts should provide MC results in HEPdata too
    # Here I just eyeballed the powheg 4F ratio
    sm = np.array(bf) * np.array([1.25, 1.09, 0.96, 1.06, 1.33])
    bf_unc = ReadDependent(vals, col=0, error=[-1], sym_errors=True)
    assert(N == len(sm) == len(bf) == len(bf_unc))

    bf = bf / sm
    bf_unc = bf_unc / sm

    cov_vals = ReadDependent(corr, col=0)
    assert((N * N) == len(cov_vals))

    cor = ROOT.TMatrixDSym(N)
    cov = ROOT.TMatrixDSym(N)

    for i in range(N):
        for j in range(N):
            cov[i][j] = cov_vals[i * N + j] / (sm[i] * sm[j])
            cor[i][j] = cov[i][j] / (bf_unc[i] * bf_unc[j])
    return {
        'N': N,
        'labels': labels,
        'sm': sm,
        'bf': bf,
        'bf_unc': bf_unc,
        'cor': cor,
        'cov': cov
    }


def MergeCov(covs=list()):
    Ntot = sum([X.GetNcols() for X in covs])
    print('Ntot = %i' % Ntot)
    cov = ROOT.TMatrixDSym(Ntot)
    pos = 0
    for c in covs:
        cov.SetSub(pos, pos, c)
        pos += c.GetNcols()
    cov.Print()
    return cov


hgg = ReadHgg()
wg = ReadWg()
singlet = ReadSingleT()

labels = hgg['labels'] + wg['labels'] + singlet['labels']
bf = np.concatenate([hgg['bf'], wg['bf'], singlet['bf']])
bf_unc = np.concatenate([hgg['bf_unc'], wg['bf_unc'], singlet['bf_unc']])

cov = MergeCov([hgg['cov'], wg['cov'], singlet['cov']])

xvars = list()
muvars = list()
xvec = ROOT.RooArgList()
mu = ROOT.RooArgList()

for l, v, u in zip(labels, bf, bf_unc):
    print(l, v, u)
    xvars.append(ROOT.RooRealVar(l, '', v, v - 10. * u, v + 10. * u))
    muvars.append(ROOT.RooRealVar(l + '_In', '', v, v - 10. * u, v + 10. * u))
    xvars[-1].setVal(1)
    muvars[-1].setConstant(True)
    xvec.add(xvars[-1])
    mu.add(muvars[-1])

xvec.Print('v')
mu.Print('v')

pdf = ROOT.RooMultiVarGaussian('pdf', '', xvec, mu, cov)
dat = ROOT.RooDataSet('global_obs', '', ROOT.RooArgSet(mu))
dat.add(ROOT.RooArgSet(mu))
nll = pdf.createNLL(dat)
minim = ROOT.RooMinimizer(nll)
minim.setEps(0.01)
# minim.setVerbose(False)
minim.minimize('Minuit2', 'migrad')
minim.hesse()
result = minim.save()
result.Print()
minim.setPrintLevel(-1)
