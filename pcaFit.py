import yaml
import json
import numpy as np
import ROOT
from array import array
from math import sqrt
import sys
import argparse
from tools import Measurement, CovTMatrix, MergeCov, ParameterUncerts

sys.path.append('EFT2Obs/scripts')
from eftscaling import EFTScaling

ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)

parser = argparse.ArgumentParser()
parser.add_argument('input', nargs='+')
args = parser.parse_args()

channel_data = dict()

# Read in the data for each channel
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

# Combine the following channels
combine_channels = channel_data.keys()

# Merge the lists of labels, best-fit values, uncertainties and covariance matrices
labels, bfs, covs = list(), list(), list()
for X in combine_channels:
    labels.extend(channel_data[X][0].bin_labels)
    bfs.extend(channel_data[X][0].bf)
    covs.append(CovTMatrix(channel_data[X][0].cov))
bf = np.array(bfs)
cov = MergeCov(covs)
bf_unc = ParameterUncerts(cov)

# Construct the model inside a workspace
w = ROOT.RooWorkspace()

# Load principal components from json file
with open('principalcomponents.json','r') as f: pca = json.load(f)

# Add a RooRealVar to the workspace for each principal component
PCs = list()
PCs_arglist = ROOT.RooArgList()
for i in range(len(pca['eigenvectors'])):
    expected_unc = 1./np.sqrt(abs(pca['eigenvalues'][i]))
    PCs.append(ROOT.RooRealVar('pc%s' % i,'pc%s' % i, -1-3*expected_unc, 1+3*expected_unc))
    PCs_arglist.add(PCs[-1])
    getattr(w, 'import')(PCs[-1])

# Define the Wilson coefficients in terms of linear combinations
# of principal components and add them to the workspace. 
# To get this system of equations, the matrix of eigenvectors has to be inverted.
eigenvectors_inv = np.linalg.inv(pca['eigenvectors'])

for i in range(len(pca['POIs'])):
    formula='%s*@0' % eigenvectors_inv[i][0]
    for j in range(1, len(eigenvectors_inv[i])):
        if eigenvectors_inv[i][j]>0: formula += '+'
        formula += '%s*@%s' % (eigenvectors_inv[i][j], j)
    getattr(w, 'import')(ROOT.RooFormulaVar(pca['POIs'][i], formula, PCs_arglist))

# Principal components with small eigenvalue (-> large uncertainty) are fixed to zero
unc_threshold = 5.
for i in range(len(pca['eigenvalues'])): 
    expected_unc = 1./np.sqrt(abs(pca['eigenvalues'][i]))
    if expected_unc > unc_threshold:
        w.var('pc%s' % i).setConstant(True)
        w.var('pc%s' % i).setVal(0.)

# List bin labels that we do not want to introduce an EFT parameterisation for
skip_bins = []

# Loop through channels in the combination, then loop through the EFT2Obs
# scaling data for each one
for X in combine_channels:
    scalings = channel_data[X][1]
    for sc in scalings:
        # Loop through scaling data (most channels only have one, but e.g. Wgamma has it split in three)
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
                # for now we only use linear terms: sum(A_i*c_i)
                if len(vars)>1: continue
                
                expr_parts.append('%g*%s' % (val, vars[0]))
            cmd = 'sum::%s(%s)' % (bin_label, ','.join(expr_parts))
            func = w.factory(cmd)


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

# Fit
minim.setEps(0.01)
minim.setVerbose(False)
minim.setPrintLevel(-1)
minim.minimize('Minuit2','migrad')

# save RooFitResult
#fitresult = minim.save()
#fitresult.Print()
#sys.exit(0)

# Keep a snapshot of initial values to be able to reset to
snapshot = {}
for i in range(len(pca['eigenvalues'])):
    snapshot['pc%s' % i] = w.var('pc%s' % i).getVal()

for i in range(len(pca['eigenvalues'])): 
    # Main loop through principal components
    expected_unc = 1./np.sqrt(abs(pca['eigenvalues'][i]))
    if expected_unc > unc_threshold: continue

    # First reset all PCs to the snapshot/best-fit value
    for j in range(len(pca['eigenvalues'])):
        w.var('pc%s' % j).setVal(snapshot['pc%s' % j])

    # The PC we're going to scan is set constant
    param = w.var('pc%s' % i)
    param.setConstant(True)
    
    # Create an output ROOT file in the same format combine makes
    # this is just so we can use existing NLL plotting script
    fout = ROOT.TFile('scan_pc%s.root' % i, 'RECREATE')
    tout = ROOT.TTree('limit', 'limit')

    a_r = array('f', [param.getVal()])
    a_deltaNLL = array('f', [0])
    a_quantileExpected = array('f', [1])

    tout.Branch('pc%s' % i, a_r, 'pc%s/f' % i)
    tout.Branch('deltaNLL', a_deltaNLL, 'deltaNLL/f')
    tout.Branch('quantileExpected', a_quantileExpected, 'quantileExpected/f')
    
    # Fill the initial (best-fit values)
    tout.Fill()
    nll0 = nll.getVal()

    # Now do a scan
    npoints = 50
    r_min = param.getMin()
    r_max = param.getMax()
    width = (float(r_max) - float(r_min)) / float(npoints)
    r = r_min + 0.5 * width

    for p in range(npoints):
        param.setVal(r)
        a_r[0] = r
        
        minim.minimize('Minuit2', 'migrad')
        nllf = nll.getVal()
        print('pc%s = %f; nll0 = %f; nll = %f, deltaNLL = %f' % (i, r, nll0, nllf,  nllf - nll0))
        a_deltaNLL[0] = nllf - nll0
        if a_deltaNLL[0] > 0.: tout.Fill()
        r += width
    tout.Write()
    fout.Close()
    param.setConstant(False)
