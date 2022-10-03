"""
Create RooWorkspaces for both the linearised and the linear+quadratic 
model and save them to a root file. Optionally a rotation matrix obtained 
by the script pca.py can be given as input to create the workspaces
in the rotated basis.

arguments:
    - input: list of input measurements and EFT scalings, e.g.
        e.g. --input ww:measurements/ATLAS_WW_parsed.yaml:scalings/ATLAS_2019_I1734263_d04-x01-y01.json \
             hgg:measurements/CMS_hgg_STXS.json:scalings/HiggsTemplateCrossSections_HTXS_stage1_2_pTjet30.json
    - rundir: directory for output (by default: rundir_{label1_label2_etc})
    - output: name of output root file (by default: workspace.root)
    - rotation: name of json file with rotation matrix (optional)
    - coeffs: list of Wilson coefficients (optional),
        will be ignored if a rotation matrix is used
        if not specified, all coefficients found in input files are used
"""
# TODO:
# - list of bins to set constant (those that are in measurement but not in scalings): make automatic

import sys
import os
import yaml
import json
import ROOT
import numpy as np
from argparse import ArgumentParser
from collections import OrderedDict
from tools import Measurement, CovTMatrix, MergeCov, ParameterUncerts
sys.path.append('EFT2Obs/scripts')
from eftscaling import EFTScaling

ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)

parser = ArgumentParser()
parser.add_argument('--input', '-i', nargs='+', 
                    help='[label]:[measurement]:[scaling]')
parser.add_argument('--rundir', '-d', default=None, 
                    help='name of run directory (output saved here)')
parser.add_argument('--output', '-o', default=None, 
                    help='name of output file')
parser.add_argument('--rotation', default=None, 
                    help='input file with rotation matrix (optional)')
parser.add_argument('--coeffs', nargs='+', default=None,
                    help='list of Wilson coefficients (optional)')
args = parser.parse_args()


# name of run directory
if args.rundir is None:
    rundir = 'rundir_%s' % '_'.join([ch.split(':')[0] for ch in args.input])
else:
    rundir = args.rundir

if not os.path.isdir(rundir):
    os.makedirs(rundir)


# name of output file
if args.output is None:
    args.output = 'workspace%s.root' % ('_warsaw' if args.rotation is None else '')
else:
    if not args.output.split('.')[-1] == 'root': 
        args.output += '.root'

outfile = os.path.join(rundir, args.output)


# look for rotation matrix input file
if args.rotation is not None:
    if not args.rotation.split('.')[-1] == 'json': 
        args.rotation += '.json'
    # if not found, check in run directory
    if not os.path.isfile(args.rotation):
        joined_path = os.path.join(rundir, args.rotation)
        if os.path.isfile(joined_path):
            args.rotation = joined_path
        else:
            sys.exit('Rotation matrix %s not found' % args.rotation)


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


# Merge the lists of bin labels, POIs, best-fit values, uncertainties and covariance matrices
cov_matrix = MergeCov([CovTMatrix(channel_data[X][0].cov) for X in channel_data])
bf = np.concatenate([channel_data[X][0].bf for X in channel_data])
bf_unc = ParameterUncerts(cov_matrix)
bin_labels = list()
POIs = list()
for X in channel_data:
    bin_labels.extend(channel_data[X][0].bin_labels)
    for scaling in channel_data[X][1]:
        for par in scaling.parameters():
            if par not in POIs:
                POIs.append(par)


# If list of Wilson coefficients is given as input, remove all other
# coefficients from the list 'POIs'. 
if args.coeffs is not None:
    if args.rotation is not None:
        print('WARNING: Using rotation matrix %s, ignore --coeffs\n' % args.rotation)
    else:
        POIs = [p for p in POIs if p in args.coeffs]
        for poi in [p for p in args.coeffs if p not in POIs]:
            print('WARNING: Wilson coefficient %s not found in input' % poi)


# Load rotation matrix from json file. If no rotation matrix is given as
# input, use identity matrix of dim {Nr. of Wilson coeffs} i.e. no rotation
if args.rotation is None:
    rot = {'xpars': [], 'ypars': POIs, 
           'rotationmatrix': np.identity(len(POIs)),
           'eigenvalues': [1. for p in POIs]}
else:
    print('>> Reading rotation matrix %s' % args.rotation)
    with open(args.rotation, 'r') as f: 
        rot = json.load(f)
print('POIs: %s\n' % ', '.join(rot['ypars']))

# Create two workspaces, one will only contain linear terms, 
# the other both linear and quadratic terms
wsp_lin = ROOT.RooWorkspace('wsp_lin')
wsp_quad = ROOT.RooWorkspace('wsp_quad')

# Add a RooRealVar to both workspaces for each rotated parameter
ypars = list()
ypars_arglist = ROOT.RooArgList()
for i,ypar in enumerate(rot['ypars']):
    expected_unc = 1./np.sqrt(abs(rot['eigenvalues'][i]))
    ypars.append(ROOT.RooRealVar(ypar, ypar, -1-3*expected_unc, 1+3*expected_unc))
    ypars_arglist.add(ypars[-1])
    getattr(wsp_lin, 'import')(ypars[-1])
    getattr(wsp_quad, 'import')(ypars[-1])


# Redefine the Wilson coefficients in terms of the rotated parameters. 
# To get this system of equations, the rotation matrix is inverted.
rotmatrix_inv = np.linalg.inv(rot['rotationmatrix'])
for i,xpar in enumerate(rot['xpars']):
    row = rotmatrix_inv[i]
    formula='%s*@0' % row[0]
    for j in range(1, len(row)):
        if row[j]>0: 
            formula += '+'
        formula += '%s*@%s' % (row[j], j)
    getattr(wsp_lin, 'import')(ROOT.RooFormulaVar(xpar, formula, ypars_arglist))
    getattr(wsp_quad, 'import')(ROOT.RooFormulaVar(xpar, formula, ypars_arglist))


# Loop through channels and add EFT parameterisations to the workspaces
for X in channel_data:
    # Loop through ET2Obs scalings of the channel
    for scaling in channel_data[X][1]:
        # Loop through bins in the scaling and construct a list
        # of terms of the form  1 + sum(A_i*c_i) + sum(B_ij*c_i*c_j)
        for ib in range(scaling.nbins):
            bin_label = scaling.bin_labels[ib]
            if bin_label not in bin_labels:  # no corresponding measurement bin
                print('Skipping EFT parameterisation for bin %s' % bin_label)
                continue
            
            # We will use the RooFit factory syntax for constructing a 
            # RooAddition as a sum of products. The parser is very fussy, 
            # and every term has to be a multiplication between exactly
            # two terms. Hence the "1*1" below.
            expr_parts_lin = ['1*1']
            expr_parts_quad = ['1*1']

            for term in scaling.terms:
                val = term.val[ib]
                vars = list(term.params)  # List of one or two coeffs
                if vars[0] not in POIs: 
                    continue

                # If two coeffs (i.e. a quadratic term), define a RooProduct for it:
                if len(vars) > 1:
                    if vars[1] not in POIs: 
                        continue
                    prodname = '_X_'.join(vars)
                    wsp_quad.factory('prod::%s(%s)' % (prodname, ','.join(vars)))
                    expr_parts_quad.append('%g*%s' % (val, prodname))
                else:
                    expr_parts_lin.append('%g*%s' % (val, vars[0]))
                    expr_parts_quad.append('%g*%s' % (val, vars[0]))
                
            wsp_lin.factory('sum::%s(%s)' % (bin_label, ','.join(expr_parts_lin)))
            wsp_quad.factory('sum::%s(%s)' % (bin_label, ','.join(expr_parts_quad)))


# xvars/x_arglist:   initially the free parameters for the fid. bins, that 
#                    will then be redefined as functions of the STXS parameters
# muvars/mu_arglist: The "global observables", representing the measured values
xvars = list()
muvars = list()
x_arglist = ROOT.RooArgList()
mu_arglist = ROOT.RooArgList()

doAsimov = False

# Right now we can only fit to a subset of the STXS bins. 
# All of these POIs will need to be set as constant.
set_constant = ["GG2H_GE2J_MJJ_0_350_PTH_0_60", "GG2H_GE2J_MJJ_0_350_PTH_60_120", "GG2H_0J_PTH_0_10", "GG2H_0J_PTH_GT10", "GG2H_1J_PTH_0_60", "GG2H_1J_PTH_120_200", "GG2H_1J_PTH_60_120", "GG2H_GE2J_MJJ_0_350_PTH_0_60   ", "GG2H_GE2J_MJJ_0_350_PTH_120_200", "GG2H_GE2J_MJJ_0_350_PTH_60_120 ", "GG2H_PTH_200_300", "GG2H_PTH_300_450", "GG2H_PTH_GT450", "QQ2HLL", "QQ2HLNU_PTV_0_75", "QQ2HLNU_PTV_75_150", "QQ2HLNU_PTV_GT150", "TH", "TTH_PTH_0_60", "TTH_PTH_120_200", "TTH_PTH_200_300", "TTH_PTH_60_120", "TTH_PTH_GT300"]

for l, v, u in zip(bin_labels, bf, bf_unc):
    if doAsimov: v = 1
    xvars.append(ROOT.RooRealVar(l, '', v, v - 10. * u, v + 10. * u))
    muvars.append(ROOT.RooRealVar(l + '_In', '', v, v - 10. * u, v + 10. * u))
    if l in set_constant: xvars[-1].setConstant(True)
    # All global observables should be set constant
    muvars[-1].setConstant(True)
    x_arglist.add(xvars[-1])
    mu_arglist.add(muvars[-1])

# Construct the pdf as a multivariate Gaussian from the covariance matrix
pdf = ROOT.RooMultiVarGaussian('pdf', '', x_arglist, mu_arglist, cov_matrix)

# Now import the pdf into the workspaces that were already populated with
# EFT scaling functions. Here we exploit "RecycleConflictNodes" - the
# functions in the wsp have the same names as the xvars that go into the PDF
# meaning they will be replaced by the existing functions.
getattr(wsp_lin, 'import')(pdf, ROOT.RooFit.RecycleConflictNodes())
getattr(wsp_quad, 'import')(pdf, ROOT.RooFit.RecycleConflictNodes())


# Save both workspaces to root file
wsp_lin.writeToFile(outfile)
wsp_quad.writeToFile(outfile, False)  # recreate=False
print('\n<< Saving workspace: %s\n' % outfile)


