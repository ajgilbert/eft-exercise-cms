import yaml
import json
import numpy as np
import ROOT
from math import sqrt
import sys
import os
import argparse
from tools import Measurement, CovTMatrix, MergeCov, ParameterUncerts
import python.plotting as plot
import python.linalg as linalg

sys.path.append('EFT2Obs/scripts')
from eftscaling import EFTScaling

ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()
plot.SetCorrMatrixPalette()

parser = argparse.ArgumentParser()
parser.add_argument('input', nargs='+', help='label:measurement:parameterisation')
parser.add_argument('--plots',  '-p', action='store_true', help='draw plots')
parser.add_argument('--output', '-o', default='plots/PCA', help='output directory for plots')
args = parser.parse_args()

if args.plots and not os.path.isdir(args.output): os.makedirs(args.output)

def TMatrixToTH2(matr, xlabels=[], ylabels=None, drawDiagonal=False):
    hist_name, hist_counter = matr.GetName(), 0 # needed to avoid memory leak when calling this function several times
    while ROOT.gROOT.FindObjectAny(hist_name+str(hist_counter)): hist_counter+=1
    
    Ncol, Mrow = matr.GetNcols(), matr.GetNrows()
    res = ROOT.TH2D(matr.GetName()+str(hist_counter), matr.GetTitle(), Ncol, 0, Ncol, Mrow, 0, Mrow)

    if ylabels is None: ylabels=xlabels[:]
    for i,xlab in enumerate(xlabels): res.GetXaxis().SetBinLabel(i+1, xlab)
    for j,ylab in enumerate(ylabels): res.GetYaxis().SetBinLabel(Mrow-j, ylab)
    
    for icol in range(1, Ncol+1):
        for jrow in range(1, Mrow+1-(icol-1) if drawDiagonal else Mrow+1):
            res.SetBinContent(icol, jrow, matr[Mrow-jrow][icol-1])
    
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

# List bin labels that we do not want to introduce an EFT parameterisation for
skip_bins = ['QQ2HQQ_FWDH', 'QQ2HQQ_0J', 'QQ2HQQ_1J', 'QQ2HQQ_GE2J_MJJ_0_60', 'QQ2HQQ_GE2J_MJJ_120_350']

# Combine the following channels
combine_channels = channel_data.keys()

# Merge the lists of labels, POIs, and covariance matrices
labels, POIs, covs = list(), list(), list()

for X in combine_channels:
    for sc in channel_data[X][1]:
        for par in sc.parameters():
            # Loops through the EFT coefficients defined for this parameterisation
            if par not in POIs:
                POIs.append(par)
    labels.extend(channel_data[X][0].bin_labels)
    covs.append(CovTMatrix(channel_data[X][0].cov))

cov = MergeCov(covs)

# matrix with linear parametrizations A_i:
# measurement bins in rows, Wilson coefficients in columns
lin_param = ROOT.TMatrixD(len(labels), len(POIs))

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
            # Now loop through bins in the EFT2Obs data and fill the linear parametrization matrix
            for term in sc.terms:
                val = term.val[ib]
                # Will be a list of one or two coeffs. If two coeffs (i.e. a quadratic term) we'll skip it
                vars = list(term.params)
                if len(vars) > 1: continue

                lin_param[labels.index(bin_label)][POIs.index(vars[0])] = val

if args.plots:
    ROOT.gStyle.SetCanvasDefW(670)
    ROOT.gStyle.SetPadRightMargin(0.12)
    ROOT.gStyle.SetPaintTextFormat('.2f')

    # plot linear parameterisation matrix
    c_lin = ROOT.TCanvas('c_lin','c_lin',1000,1000)
    h_lin = TMatrixToTH2(lin_param, POIs, labels)
    h_lin.SetMarkerSize(0.5)
    h_lin.GetXaxis().SetLabelSize(0.02)
    h_lin.GetYaxis().SetLabelSize(0.02)
    h_lin.Draw('colz text')
    c_lin.Print(os.path.join(args.output, 'linparam_matrix.pdf'))

# Fisher information matrix
fisher = linalg.TMMultiply(linalg.TMTranspose(lin_param), linalg.TMInvert(cov), lin_param)

if args.plots:
    # plot fisher matrix
    c_fisher = ROOT.TCanvas('c_fisher','c_fisher',1000,1000)
    h_fisher = TMatrixToTH2(fisher, POIs)
    h_fisher.SetMarkerSize(0.5)
    h_fisher.GetXaxis().SetLabelSize(0.025)
    h_fisher.GetYaxis().SetLabelSize(0.025)
    h_fisher.Draw('colz text')
    c_fisher.Print(os.path.join(args.output, 'fishermatrix.pdf'))

    # plot normalized fisher matrix (A_ij -> A_ij/sqrt(A_ii*A_jj))
    c_fishernorm = ROOT.TCanvas('c_fishernorm','c_fishernorm',1000,1000)
    h_fishernorm = TMatrixToTH2(linalg.CovToCorr(fisher), POIs)
    h_fishernorm.SetMarkerSize(0.5)
    h_fishernorm.GetXaxis().SetLabelSize(0.025)
    h_fishernorm.GetYaxis().SetLabelSize(0.025)
    h_fishernorm.Draw('colz text')
    c_fishernorm.Print(os.path.join(args.output, 'fishermatrix_norm.pdf'))

# get eigenvalues and eigenvectors of Fisher matrix
eigen_fisher = ROOT.TMatrixDEigen(fisher)
eigenvalues, eigenvectors = eigen_fisher.GetEigenValues(), eigen_fisher.GetEigenVectors()

if args.plots:
    # plot matrix with eigenvectors of Fisher matrix (principal components) in columns
    c_eigvectors = ROOT.TCanvas('c_eigvectors','c_eigvectors',1000,1000)
    h_eigvectors = TMatrixToTH2(eigenvectors, ['v_{%s}' % i for i in range(eigenvalues.GetNcols())], POIs)
    h_eigvectors.SetMarkerSize(0.5)
    h_eigvectors.Draw('colz text')
    c_eigvectors.Print(os.path.join(args.output, 'eigenvectors.pdf'))

    # plot matrix with principal components in rows, but only those with
    # expected uncertainty smaller than some threshold
    ylabels, PCs = [], []
    unc_threshold = 5.

    for i in range(eigenvalues.GetNcols()):
        # expected uncertainty of measurement in direction of an eigenvector ~ 1/sqrt(eigenvalue)
        unc = 1./np.sqrt(abs(eigenvalues[i][i]))
        if unc < unc_threshold:
            PCs.append(i)
            ylabels.append('PC%s (\sigma = %.3f)' % (i, unc)) 

    h_pca = TMatrixToTH2(linalg.TMSubmatrix( linalg.TMTranspose(eigenvectors), rows=PCs ), POIs, ylabels)
    rpc = float(h_pca.GetNbinsY())/h_pca.GetNbinsX()
    c_pca = ROOT.TCanvas('c_pca','c_pca',1000,int(1000*rpc))
    h_pca.SetMarkerSize(0.5/rpc)
    h_pca.GetXaxis().SetLabelSize(0.02/rpc)
    h_pca.GetYaxis().SetLabelSize(0.02/rpc)
    h_pca.Draw('colz text')
    c_pca.Print(os.path.join(args.output, 'pca.pdf'))

# save to json
with open('principalcomponents.json', 'w') as f:
    json.dump({'POIs': POIs, 
               'eigenvalues': [eigenvalues[i][i] for i in range(eigenvalues.GetNrows())], 
               'eigenvectors': [list(linalg.TMatrixToArray(eigenvectors).T[i]) for i in range(eigenvectors.GetNcols())]}, 
               f, indent=2)


