"""
Principal Component Analysis (PCA)

- Get Fisher information matrix in EFT basis: 
        F = P^{T} H P
    - H: block-diagonal Hessian matrix of the measurements
    - P: linear parameterisation matrix

- Eigendecomposition of the Fisher information matrix
    - Eigenvectors x_i of F are othogonal directions in Wilson coefficient space
    - Eigenvalues lambda_i related to uncertainty of measurement in direction of 
        eigenvector x_i:  sigma = 1/sqrt(lambda)

- Save Fisher matrix, rotation matrix (x_1 ... x_N)^T, and eigenvalues to json 

arguments:
    - input: list of input measurements and EFT scalings, e.g.
        e.g. --input ww:measurements/ATLAS_WW_parsed.yaml:scalings/ATLAS_2019_I1734263_d04-x01-y01.json \
             hgg:measurements/CMS_hgg_STXS.json:scalings/HiggsTemplateCrossSections_HTXS_stage1_2_pTjet30.json
    - rundir: directory for output (by default: rundir_{label1_label2_etc})
    - output: name of output json file (by default: rotmatrix.json)
"""

import sys
import os
import yaml
import json
import ROOT
import numpy as np
from argparse import ArgumentParser
from collections import OrderedDict
from tools import Measurement, CovTMatrix, MergeCov
from python.linalg import TMMultiply, TMTranspose, TMInvert, TMatrixToArray
sys.path.append('EFT2Obs/scripts')
from eftscaling import EFTScaling


parser = ArgumentParser()
parser.add_argument('--input', '-i', nargs='+', 
                    help='[label]:[measurement]:[scaling]')
parser.add_argument('--rundir', '-d', default=None, 
                    help='name of run directory (output saved here)')
parser.add_argument('--output', '-o', default='rotmatrix.json', 
                    help='name of output file')
args = parser.parse_args()


# name of run directory
if args.rundir is None:
    rundir = 'rundir_%s' % '_'.join([ch.split(':')[0] for ch in args.input])
else:
    rundir = args.rundir

if not os.path.isdir(rundir):
    os.makedirs(rundir)


# name of output file 
if not args.output.split('.')[-1] == 'json': 
    args.output += '.json'

outfile = os.path.join(rundir, args.output)


# Read in the data for each channel
print('')
channel_data = OrderedDict()
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
print('')


# Merge the lists of bin labels, POIs
bin_labels = list()
POIs = list()
for X in channel_data:
    bin_labels.extend(channel_data[X][0].bin_labels)
    for scaling in channel_data[X][1]:
        for par in scaling.parameters():
            if par not in POIs:
                POIs.append(par)

# construct block-diagonal covariance matrix
cov_matrix = TMatrixToArray(MergeCov([CovTMatrix(channel_data[X][0].cov) for X in channel_data]))

# make sure the combined covariance matrix is not rank-deficient (-> can be inverted)
assert(cov_matrix.shape[0] == cov_matrix.shape[1] == np.linalg.matrix_rank(cov_matrix))

# matrix with linear parametrisations A_i:
# measurement bins in rows, Wilson coefficients in columns
lin_param = np.zeros([len(bin_labels), len(POIs)])

# Loop through channels and fill the linear parametrisation matrix
for X in channel_data:
    # Loop through ET2Obs scalings of the channel
    for scaling in channel_data[X][1]:
        # Loop through bins in the scaling
        for ib in range(scaling.nbins):
            bin_label = scaling.bin_labels[ib]
            if bin_label not in bin_labels:  # no corresponding measurement bin
                print('Skipping EFT parameterisation for bin %s' % bin_label)
                continue
            
            for term in scaling.terms:
                val = term.val[ib]
                vars = list(term.params)  # List of one or two coeffs
                if len(vars) == 1:  # linear term
                    lin_param[bin_labels.index(bin_label)][POIs.index(vars[0])] = val


# Rotate the inverted covariance matrix to the EFT basis
fisher_matrix = np.linalg.multi_dot((lin_param.T, np.linalg.pinv(cov_matrix), lin_param))

# Eigendecomposition of the Fisher matrix
eigen_fisher = np.linalg.svd(fisher_matrix)
eigenvalues = eigen_fisher[1]
eigenvectors = eigen_fisher[0].T

# Save matrices to json file
print('\n<< Saving matrices: %s\n' % outfile)
with open(outfile, 'w') as fout:
    json.dump(OrderedDict([
        ('input', args.input),
        ('xpars', POIs), 
        ('ypars', ['EV%s' % (i+1) for i in range(len(eigenvalues))]),  # ypars = rotationmatrix * xpars
        ('fishermatrix', fisher_matrix.tolist()),
        ('rotationmatrix', eigenvectors.tolist()),
        ('covmatrix', cov_matrix.tolist()),
        ('eigenvalues', eigenvalues.tolist())
        ]), fout, indent=2)


