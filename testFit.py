import yaml
import json
import numpy as np
import ROOT
from array import array
from math import sqrt
import sys
import argparse
from tools import Measurement

sys.path.append('EFT2Obs/scripts')
from eftscaling import EFTScaling

ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)

parser = argparse.ArgumentParser()
parser.add_argument('input', nargs='+')
args = parser.parse_args()


def CovTMatrix(cov):
    shape = np.shape(cov)
    assert(shape[0] == shape[1])
    N = shape[0]
    res = ROOT.TMatrixDSym(N)
    for i in range(N):
        for j in range(N):
            res[i][j] = cov[i][j]
    return res

def MergeCov(covs=list()):
    # Input: list of TMatrixD
    # Output: block-diagonal TMatrixD from combination of inputs
    Ntot = sum([X.GetNcols() for X in covs])
    cov = ROOT.TMatrixDSym(Ntot)
    pos = 0
    for c in covs:
        cov.SetSub(pos, pos, c)
        pos += c.GetNcols()
    # cov.Print()
    return cov

def ParameterUncerts(matrix):
    N = matrix.GetNcols()
    res = []
    for i in range(N):
        res.append(sqrt(matrix[i][i]))
    return res

channel_data = dict()

for channel in args.input:
    label, measurement_file, param_files = channel.split(':')
    print('>> Reading measurement {} and parameterisation {}'.format(measurement_file, param_files))
    if measurement_file.split('.')[-1]=='yaml':
      measurement = Measurement.fromYAML(measurement_file)
    elif measurement_file.split('.')[-1]=='json':
      measurement = Measurement.fromJSON(measurement_file)
    else:
      print('File format not supported')
    params = [EFTScaling.fromJSON(X) for X in param_files.split(',')]
    channel_data[label] = (measurement, params)



# Read in the data for each channel

# Combine the following channels
combine_channels = channel_data.keys()

# Merge the lists of labels, best-fit values, uncertainties and covariance matrices
labels = list()
for X in combine_channels:
    labels.extend(channel_data[X][0].bin_labels)
bf = np.concatenate([channel_data[X][0].bf for X in combine_channels])
cov = MergeCov([CovTMatrix(channel_data[X][0].cov) for X in combine_channels])
bf_unc = ParameterUncerts(cov)
# cov.Print()

# sys.exit(0)
# Construct the model inside a workspace
w = ROOT.RooWorkspace()

# TODO: for now we hardcode reasonable ranges for each POI, but this
# should be set via some external config file
scan_ranges = {
    "chdd": [-5, 5],
    "chj3": [-1, 1],
    "chl3": [-3, 3],
    "chwb": [-1, 1],
    "clj3": [-0.4, 0.4],
    "cll1": [-3, 3],
    "cw": [-0.1, 0.1],
    "chg": [-1, 1],
    "chb": [-5, 5],
    "chbox": [-10, 10],
    "chd": [-5, 5],
    "chj1": [-5, 5],
    "chu": [-5, 5],
    "chw": [-3, 3],
    "cbgre": [-30, 30],
    "cbwre": [-10, 10],
    "chq3": [-4, 4],
    "chtbre": [-30, 30],
    "cqj31": [-1, 1],
    "cqj38": [-1, 1],
    "ctgre": [-10, 10],
    "ctwre": [-3, 3]
}

POIs = list()

# List bin labels that we do not want to introduce an EFT parameterisation for
skip_bins = []


def makePOI(name, range_dict, wsp, default_range=[-1., 1.]):
    poi_range = default_range if name not in range_dict else range_dict[name]
    wsp.factory('%s[0,%g,%g]' % (name, poi_range[0], poi_range[1]))




# Loop through channels in the combination, then loop through the EFT2Obs
# scaling data for each one
for X in combine_channels:
    scalings = channel_data[X][1]
    for sc in scalings:
        # Loop through scaling data (most channels only have one, but e.g. Wgamma has it split in three)
        for par in sc.parameters():
            # Loops through the EFT coefficients defined for this parameterisation
            if par not in POIs:
                POIs.append(par)
                makePOI(par, scan_ranges, w)
        for ib in range(sc.nbins):
            bin_label = sc.bin_labels[ib]
            if bin_label in skip_bins:
                print('Skipping EFT parameterisation for bin %s' % bin_label)
                continue
            # Now loop through bins in the EFT2Obs data, construct a list of terms of the form:
            # 1 + sum(A_i*c_i) + sum(B_ij*c_i*c_j)
            # We will use the RooFit factory syntax for constructing a RooAddition as a sum of 
            # products. The parser is very fussy, and every term has to be a multiplication between
            # exactly two terms. Hence the "1*1" below.
            expr_parts = ['1*1']
            for term in sc.terms:
                val = term.val[ib]
                # Will be a list of one or two coeffs
                vars = list(term.params)

                # If two coeffs (i.e. a quadratic term) we'll define a RooProduct for it:
                if len(vars) > 1:
                    prodname = '_X_'.join(vars)
                    w.factory('prod::%s(%s)' % (prodname, ','.join(vars)))
                    expr_parts.append('%g*%s' % (val, prodname))
                else:
                    expr_parts.append('%g*%s' % (val, vars[0]))
                # print(val, vars)
            cmd = 'sum::%s(%s)' % (bin_label, ','.join(expr_parts))
            # print(cmd)
            func = w.factory(cmd)
            # func.Print()
        

# xvars/xvec: initially the free parameters for the fid. bins, that 
#             will then be redefined as functions of the STXS parameters
# muvars/mu:  The "global observables", representing the measured values
xvars = list()
muvars = list()
xvec = ROOT.RooArgList()
mu = ROOT.RooArgList()

doAsimov = False

# Right now we can only fit to a subset of the STXS bins. All of these POIs will need to be set as constant.
set_constant = ["GG2H_GE2J_MJJ_0_350_PTH_0_60", "GG2H_GE2J_MJJ_0_350_PTH_60_120", "GG2H_0J_PTH_0_10", "GG2H_0J_PTH_GT10", "GG2H_1J_PTH_0_60", "GG2H_1J_PTH_120_200", "GG2H_1J_PTH_60_120", "GG2H_GE2J_MJJ_0_350_PTH_0_60   ", "GG2H_GE2J_MJJ_0_350_PTH_120_200", "GG2H_GE2J_MJJ_0_350_PTH_60_120 ", "GG2H_PTH_200_300", "GG2H_PTH_300_450", "GG2H_PTH_GT450", "QQ2HLL", "QQ2HLNU_PTV_0_75", "QQ2HLNU_PTV_75_150", "QQ2HLNU_PTV_GT150", "TH", "TTH_PTH_0_60", "TTH_PTH_120_200", "TTH_PTH_200_300", "TTH_PTH_60_120", "TTH_PTH_GT300"]

for l, v, u in zip(labels, bf, bf_unc):
    # print(l, v, u)
    if doAsimov:
        v = 1
    xvars.append(ROOT.RooRealVar(l, '', v, v - 10. * u, v + 10. * u))
    muvars.append(ROOT.RooRealVar(l + '_In', '', v, v - 10. * u, v + 10. * u))
    if l in set_constant:
        xvars[-1].setConstant(True)
    # All global observables should be set constant
    muvars[-1].setConstant(True)
    xvec.add(xvars[-1])
    mu.add(muvars[-1])

xvec.Print('v')
mu.Print('v')

# Construct the pdf as a multivariate Gaussian from the covariance matrix
pdf = ROOT.RooMultiVarGaussian('pdf', '', xvec, mu, cov)
# Now import the pdf into the workspace that was already populated with
# EFT scaling functions. Here we exploit "RecycleConflictNodes" - the
# functions in the wsp have the same names as the xvars that go into the PDF
# meaning they will be replaced by the existing functions
getattr(w, 'import')(pdf, ROOT.RooFit.RecycleConflictNodes())

# Take the re-parameterised pdf from the workspace
pdf = w.pdf('pdf')

# Define the data using the global observables
dat = ROOT.RooDataSet('global_obs', '', ROOT.RooArgSet(mu))
dat.add(ROOT.RooArgSet(mu))

# Create the NLL function and the minimizer
nll = pdf.createNLL(dat)
minim = ROOT.RooMinimizer(nll)

# Since fit is trivial, set a tight tolerance
minim.setEps(0.01)
minim.setVerbose(False)
minim.setPrintLevel(-1)

# Keep a snapshot of initial values to be able to reset to
snapshot = {}
for POI in POIs:
    snapshot[POI] = w.var(POI).getVal()
    # Start by setting all POIs constant
    w.var(POI).setConstant(True)

for POI in POIs:
    # Main loop through POIs
    # First set all POIs constant and to the snapshot value
    for POI2 in POIs:
        w.var(POI2).setConstant(True)
        w.var(POI2).setVal(snapshot[POI])

    # The POI we're going to fit:
    param = w.var(POI)
    # Float just this one, and do the fit
    param.setConstant(False)
    minim.minimize('Minuit2', 'migrad')

    # Create an output ROOT file in the same format combine
    # makes - this is just so we can use existing NLL plotting
    # script
    fout = ROOT.TFile('scan_%s.root' % (POI), 'RECREATE')
    tout = ROOT.TTree('limit', 'limit')

    a_r = array('f', [param.getVal()])
    a_deltaNLL = array('f', [0])
    a_quantileExpected = array('f', [1])

    tout.Branch(POI, a_r, '%s/f' % POI)
    tout.Branch('deltaNLL', a_deltaNLL, 'deltaNLL/f')
    tout.Branch('quantileExpected', a_quantileExpected, 'quantileExpected/f')
    # Fill the initial (best-fit values)
    tout.Fill()
    nll0 = nll.getVal()

    # Now do a scan
    npoints = 50
    param.setConstant(True)
    r_min = scan_ranges[POI][0]
    r_max = scan_ranges[POI][1]
    width = (float(r_max) - float(r_min)) / float(npoints)
    r = r_min + 0.5 * width
    # print(r, width)
    for p in range(npoints):
        param.setVal(r)
        a_r[0] = r
        # Won't be any actual fitting to do here (no floating params), but later we may float more...
        minim.minimize('Minuit2', 'migrad')
        nllf = nll.getVal()
        print('%s = %f; nll0 = %f; nll = %f, deltaNLL = %f' % (POI, r, nll0, nllf,  nllf - nll0))
        a_deltaNLL[0] = nllf - nll0
        if a_deltaNLL[0] > 0.: tout.Fill()
        r += width
    tout.Write()
    fout.Close()
