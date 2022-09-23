"""
Take workspace from root file, fit, and save fit result to json file. With the
option --scan or --scan_fixed, nll scans for each parameter are saved to root
files. Parameters with large uncertainty are fixed to zero (by default > 10).

arguments:
    - input: name of input root file containing the workspace
    - rundir: directory for input/output (by default: rundir_{N})
    - output: name of output json file (by default: fitresult_{lin/quad}.json)
    - quadratic: use linear+quadratic model (by default: linearised)
    - fix: list of parameters that should be fixed to zero (optional)
    - max_unc: parameters with uncertainty > max_unc are fixed to zero
    - scan: do nll scans for each parameter with all other parameters floating
    - scan_fixed: do nll scans for each parameter with all other parameters fixed to zero
    - print_scan: print results of parameter scans
"""

import sys
import os
import json
import ROOT
import numpy as np
from argparse import ArgumentParser
from collections import OrderedDict
from array import array
from python.linalg import TMatrixToArray

ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)

parser = ArgumentParser()
parser.add_argument('--input', '-i', default=None,
                    help='input root file with RooWorkspace')
parser.add_argument('--rundir', '-d', default=None, 
                    help='name of run directory (output saved here)')
parser.add_argument('--output', '-o', default=None, 
                    help='name of output file for fit result')
parser.add_argument('--quadratic', '-q', action='store_true', 
                    help='fit linear and quadratic terms (default: linear only)')
parser.add_argument('--max_unc', default=1./np.sqrt(0.01), 
                    help='fix all POIs with uncertainty above this value to zero')
parser.add_argument('--fix', nargs='+', default=[], 
                    help='list of POIs that should be fixed to zero (optional)')
parser.add_argument('--scan', action='store_true', 
                    help='do nll scans for each parameter')
parser.add_argument('--scan_fixed', action='store_true', 
                    help='do nll scans for each parameter with all other parameters fixed to zero')
parser.add_argument('--scan_only', action='store_true', 
                    help='skip simultaneous fit and only do nll scans')
parser.add_argument('--print_scan', action='store_true', 
                    help='print results of parameter scans')
args = parser.parse_args()

if args.input is None:
    sys.exit('no input file given (--input NAME_OF_FILE.root)')


# name of run directory
if args.rundir is None:
    if os.path.split(args.input)[0] != '':
        args.rundir = os.path.split(args.input)[0]
    else:
        i_dir = 1
        while 'rundir_%s' % i_dir in os.listdir('.'):
            i+=1
        args.rundir = 'rundir_%s' % i_dir


# look for input file
if not args.input.split('.')[-1] == 'root': 
    args.input += '.root'
# if not found, check in run directory
if not os.path.isfile(args.input):
    joined_path = os.path.join(args.rundir, args.input)
    if os.path.isfile(joined_path):
        args.input = joined_path
    else:
        sys.exit('Input file %s not found' % args.input)


# name of output file
if args.output is None:
    outfile = 'fitresult%s_%s.json' % (
        '_warsaw' if 'warsaw' in args.input else '',
        'quad' if args.quadratic else 'lin')
else:
    outfile = args.output    
    if not outfile.split('.')[-1] == 'json': 
        outfile += '.json'

outfile = os.path.join(args.rundir, outfile)


# make sure output directories exist
if not os.path.isdir(args.rundir):
    os.makedirs(args.rundir)

if args.scan or args.scan_fixed or args.scan_only:
    scandir = os.path.join(args.rundir, 'scans')
    if not os.path.isdir(scandir):
        os.makedirs(scandir)


# get workspaces from input root file
wsp_name = 'wsp_quad' if args.quadratic else 'wsp_lin'
fin = ROOT.TFile(args.input, 'READ')
wsp = fin.Get(wsp_name)
print('\n>> Taking workspace %s from %s\n' % (wsp_name, args.input))

# Take the pdf from the workspace
pdf = wsp.pdf('pdf')

# Get list of variables in the workspace
varnames = wsp.allVars().contentsString().split(',')

# Take set of global observables from the workspace
mu_argset = ROOT.RooArgSet()
poi_names = list()
for name in varnames:
    if '_In' in name:
        mu_argset.add(wsp.var(name))

# Define the data using the global observables
dat = ROOT.RooDataSet('global_obs', '', ROOT.RooArgSet(mu_argset))
dat.add(ROOT.RooArgSet(mu_argset))

# Take set of parameters from the workspace
poi_names = [v for v in varnames if not wsp.var(v).getAttribute('Constant')]

# Parameters specified by --fix are fixed to zero
for poi in poi_names: 
    if poi in args.fix:
        wsp.var(poi).setVal(0.)
        wsp.var(poi).setConstant(True)
if len(args.fix) > 0:
    print('These POIs are fixed to zero: %s\n' % ', '.join(args.fix))
poi_names = [poi for poi in poi_names if not poi in args.fix]

# Create the NLL function and the minimizer
nll = pdf.createNLL(dat)
minim = ROOT.RooMinimizer(nll)

# Fit
minim.setEps(0.01)
minim.setVerbose(False)
minim.setPrintLevel(-1)

if not args.scan_only:
    minim.minimize('Minuit2','migrad')

    # Parameters with uncertainty > args.max_unc are fixed to zero
    fixzero = list()
    for poi in poi_names:
        err = wsp.var(poi).getError()
        if err > args.max_unc:
            wsp.var(poi).setVal(0.)
            wsp.var(poi).setConstant(True)
            fixzero.append(poi)
    poi_names = [poi for poi in poi_names if not poi in fixzero]

    # Repeat fit with reduced number of parameters
    if len(fixzero) > 0:
        minim.minimize('Minuit2','migrad')
        print('\nThese POIs are fixed to zero because their uncertainty is > %.2f: %s' % 
            (args.max_unc, ', '.join(fixzero)))
        print('The maximum uncertainty can be set with the option --max_unc')

    # save RooFitResult
    fitresult = minim.save()
    fitresult.Print()


# Keep a snapshot of initial values to be able to reset to
snapshot = OrderedDict()
for poi in poi_names:
    snapshot[poi] = [wsp.var(poi).getVal(), wsp.var(poi).getError()]


if not args.scan_only:
    # Save fit result and correlation matrix to json file
    with open(outfile, 'w') as f:
        json.dump(OrderedDict([
            ('workspace', '%s: %s' % (args.input, wsp_name)), 
            ('fitresult', snapshot),
            ('correlationmatrix', TMatrixToArray(fitresult.correlationMatrix()).tolist())
            ]), f, indent=2)
    print('\n<< saving fit results to %s' % outfile)


# Likelihood scan for each parameter
if not (args.scan or args.scan_fixed or args.scan_only):
    print('')
    sys.exit(0)

outfiles = list()
for poi in poi_names: 
    # First reset all parameters to the snapshot/best-fit value
    for poi2 in poi_names:
        if args.scan_fixed:
            wsp.var(poi2).setVal(0.)
            wsp.var(poi2).setConstant(True)
        else:
            wsp.var(poi2).setVal(snapshot[poi2][0])
            wsp.var(poi2).setError(snapshot[poi2][1])

    # The parameter we are going to scan
    param = wsp.var(poi)

    # find best-fit value when all other parameters fixed to zero
    if args.scan_fixed or args.scan_only:
        param.setConstant(False)
        minim.minimize('Minuit2', 'migrad')

    param.setConstant(True)
    if args.print_scan: 
        print('\nScanning %s = %.3f +/- %.3f, all other POIs %s\n' % (poi, 
            param.getVal(), param.getError(), 
            'fixed to zero' if args.scan_fixed else 'floating'))

    # Create an output ROOT file in the same format combine makes
    # this is just so we can use existing NLL plotting script
    fname = 'scan_%s_%s' % (poi, 'quad' if args.quadratic else 'lin')
    if args.scan_fixed or args.scan_only: 
        fname += '_fixed'
    fname += '.root'
    outfiles.append(os.path.join(scandir, fname))

    fout = ROOT.TFile(os.path.join(scandir, fname), 'RECREATE')
    tout = ROOT.TTree('limit', 'limit')

    a_r = array('f', [param.getVal()])
    a_deltaNLL = array('f', [0])
    a_quantileExpected = array('f', [1])

    tout.Branch(poi, a_r, '%s/f' % poi)
    tout.Branch('deltaNLL', a_deltaNLL, 'deltaNLL/f')
    tout.Branch('quantileExpected', a_quantileExpected, 'quantileExpected/f')
    
    # Fill the initial (best-fit values)
    tout.Fill()
    nll0 = nll.getVal()

    # Now do a scan
    npoints = 50
    r_min = param.getVal()-4*param.getError()
    r_max = param.getVal()+4*param.getError()
    width = (float(r_max) - float(r_min)) / float(npoints)
    r = r_min + 0.5 * width

    for p in range(1,npoints+1):
        param.setVal(r)
        a_r[0] = r
        minim.minimize('Minuit2', 'migrad')
        nllf = nll.getVal()
        if args.print_scan: 
            print('%s: %s = %f; nll0 = %f; nll = %f, deltaNLL = %f' % (p, poi, r, nll0, nllf,  nllf - nll0))
        a_deltaNLL[0] = nllf - nll0
        if a_deltaNLL[0] > 0.: 
            tout.Fill()
        r += width
    
    tout.Write()
    fout.Close()

    param.setConstant(False)

print('')
for outf in outfiles:
    print('<< saving parameter scan %s' % outf)
print('')


