import yaml
import json
import numpy as np
import ROOT
from array import array
from math import sqrt

ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)

def ReadIndependent(entry, col=0):
    # Extract the bin labels / bin ranges from the hepData YAML
    values = entry['independent_variables'][col]['values']
    if 'value' in values[0]:
        return [X['value'] for X in values]
    else:
        return [(X['low'], X['high']) for X in values]


def ReadDependent(entry, col=0, error=list(), sym_errors=True):
    # Extract the central values or the uncertainties from a column in the
    # hepData YAML. To extract the errors, supply a list of error indicies
    # that should be summed in quadrature. Errors are then symmeterised by
    # default.
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
    # Read the HIG-19-015 "minimal merging" STXS results and correlation matrix
    with open("hepdata_inputs/HEPData-ins1851456-v2-STXS_stage_1.2_minimal_merging_scheme.yaml", "r") as f:
        vals = yaml.safe_load(f)
    with open("hepdata_inputs/HEPData-ins1851456-v2-Correlations__STXS_stage_1.2_minimal_merging_scheme.yaml", "r") as f:
        corr = yaml.safe_load(f)

    # Read the labels from the hepData table. However, for now we'll override these manually
    # with labels that correspond to the bin naming in the STXS Rivet ouput from EFT2Obs.
    labels = ReadIndependent(vals, col=0)
    # Some labels combine multiple STXS bins, that are merged in the analysis. These are given
    # as comments below.
    new_labels = [
        "GG2H_0J_PTH_0_10",
        "GG2H_0J_PTH_GT10",
        "GG2H_1J_PTH_0_60",
        "GG2H_1J_PTH_60_120",
        "GG2H_1J_PTH_120_200",
        "GG2H_GE2J_MJJ_0_350_PTH_0_60",
        "GG2H_GE2J_MJJ_0_350_PTH_60_120",
        "GG2H_GE2J_MJJ_0_350_PTH_120_200",
        "GG2H_PTH_200_300",
        "GG2H_PTH_300_450",
        "GG2H_PTH_GT450", # "GG2H_PTH_450_650", "GG2H_PTH_GT650",
        "QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25", # + "GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25",
        "QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25", # + "GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25",
        "QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25", # + "GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25",
        "QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25", # + "GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25",
        "QQ2HQQ_GE2J_MJJ_60_120",
        "QQ2HQQ_GE2J_MJJ_GT350_PTH_GT200",
        "QQ2HLNU_PTV_0_75",
        "QQ2HLNU_PTV_75_150",
        "QQ2HLNU_PTV_GT150", #"QQ2HLNU_PTV_150_250_0J",    "QQ2HLNU_PTV_150_250_GE1J", "QQ2HLNU_PTV_GT250", 
        "QQ2HLL", # "QQ2HLL_PTV_0_75", "QQ2HLL_PTV_75_150", "QQ2HLL_PTV_150_250_0J", "QQ2HLL_PTV_150_250_GE1J", "QQ2HLL_PTV_GT250", "GG2HLL_FWDH", "GG2HLL_PTV_0_75", "GG2HLL_PTV_75_150", "GG2HLL_PTV_150_250_0J", "GG2HLL_PTV_150_250_GE1J", "GG2HLL_PTV_GT250",
        "TTH_PTH_0_60",
        "TTH_PTH_60_120",
        "TTH_PTH_120_200",
        "TTH_PTH_200_300",
        "TTH_PTH_GT300",
        "TH"
    ]
    N = len(labels)

    # Just check that new_labels is the same length as the labels from the hepData,
    # before replacing them
    assert(N == len(new_labels))
    labels = new_labels

    # Read the SM cross sections from the first column
    sm = ReadDependent(vals, col=0)
    # Read the best-fit values and uncertainties from the second column
    bf = ReadDependent(vals, col=1)
    # The total errors are already computed
    bf_unc = ReadDependent(vals, col=1, error=[0], sym_errors=True)
    assert(N == len(sm) == len(bf) == len(bf_unc))

    # The above values are absolute cross sections, it's easiest for us to deal
    # with everything relative to the SM values:
    bf = bf / sm
    bf_unc = bf_unc / sm

    # Read the correlation matrix
    corr_vals = ReadDependent(corr, col=0)
    assert((N * N) == len(corr_vals))

    cor = ROOT.TMatrixDSym(N)
    cov = ROOT.TMatrixDSym(N)

    # This part is analysis dependent - have to know the (i,j) index ordering
    # in the 1D array of values
    for i in range(N):
        for j in range(N):
            cor[i][j] = corr_vals[i * N + j]
            # compute the covariance
            cov[i][j] = corr_vals[i * N + j] * bf_unc[i] * bf_unc[j]
    
    # Load the EFT2Obs scaling json
    with open('eft2obs_inputs/HiggsTemplateCrossSections_HTXS_stage1_2_pTjet30.json') as jsonfile:
        scaling = json.load(jsonfile)
    return {
        'N': N,
        'labels': labels,
        'sm': sm,
        'bf': bf,
        'bf_unc': bf_unc,
        'cor': cor,
        'cov': cov,
        'scaling': [scaling]
    }


def ReadWg():
    # Read the SMP-20-005 pTgamma x |phi_f| cross section

    # The unrolling of the 2D is phi first, then pT. To match up with the cov.
    # matrix it will be easiest to covert to pT,phi unrolled instead.
    def Reshape(arr):
        # BAD: we hardcode "3" as the number of phi bins
        return arr.reshape((3, int(len(arr) / 3))).T.reshape((-1))

    with open("hepdata_inputs/HEPData-ins1978840-v1-p_{mathrm{T}}^{gamma}_times_|phi_{f}|_cross_section.yaml", "r") as f:
        vals = yaml.safe_load(f)
    with open("hepdata_inputs/HEPData-ins1978840-v1-p_{mathrm{T}}^{gamma}_times_|phi_{f}|_correlations.yaml", "r") as f:
        corr = yaml.safe_load(f)

    labels_i = ReadIndependent(vals, col=1)
    labels = list()
    N = len(labels_i)
    for i in range(int(N / 3)):
        for j in range(3):
            labels.append('WG_pt_%i_phi_%i' % (i, j))
    sm = Reshape(ReadDependent(vals, col=1))
    bf = Reshape(ReadDependent(vals, col=0))
    bf_unc = Reshape(ReadDependent(vals, col=0, error=[0, 1], sym_errors=True))
    assert(N == len(sm) == len(bf) == len(bf_unc))

    bf = bf / sm
    bf_unc = bf_unc / sm

    corr_vals = ReadDependent(corr, col=0)

    cor = ROOT.TMatrixDSym(N)
    cov = ROOT.TMatrixDSym(N)

    # Ordering of terms in the 1D correlation array is non-intuitive,
    # compared to Hgg  above, hence the complicated index logic:
    idx = 0
    for i in range(N):
        for j in range(N - i):
            cor[i][N - j - 1] = corr_vals[idx]
            cor[N - j - 1][i] = cor[i][N - j - 1] 
            cov[i][N - j - 1] = corr_vals[idx] * bf_unc[i] * bf_unc[N - j - 1]
            cov[N - j - 1][i] = cov[i][N - j - 1]
            idx += 1

    # Rivet routine gives us three separate 1D histograms - one for each phi bin
    scalings = []
    for file in ['CMS_2021_PAS_SMP_20_005_d54-x01-y01.json', 'CMS_2021_PAS_SMP_20_005_d55-x01-y01.json', 'CMS_2021_PAS_SMP_20_005_d56-x01-y01.json']:
        with open('eft2obs_inputs/%s' % file) as jsonfile:
            scalings.append(json.load(jsonfile))

    return {
        'N': N,
        'labels': labels,
        'sm': sm,
        'bf': bf,
        'bf_unc': bf_unc,
        'cor': cor,
        'cov': cov,
        'scaling': scalings
    }


def ReadSingleT():
    with open("hepdata_inputs/HEPData-ins1744604-v1-Table_1.yaml", "r") as f:
        vals = yaml.safe_load(f)
    with open("hepdata_inputs/HEPData-ins1744604-v1-Table_2.yaml", "r") as f:
        corr = yaml.safe_load(f)

    labels = ReadIndependent(vals, col=0)
    labels = ['pt_t_bin_%i' % X for X in range(len(labels))]
    N = len(labels)
    bf = ReadDependent(vals, col=0)
    # TODO: Not able to find numerical values of MC predictions
    # Here I just eyeballed the powheg 4F ratio:
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
            # Here we are given the covariance of the differential cross section
            # directly - so we need to convert to ratio wrt. SM
            cov[i][j] = cov_vals[i * N + j] / (sm[i] * sm[j])
            cor[i][j] = cov[i][j] / (bf_unc[i] * bf_unc[j])
    
    with open('eft2obs_inputs/CMS_2019_I1744604_d13-x01-y01.json') as jsonfile:
        scaling = json.load(jsonfile)
    
    return {
        'N': N,
        'labels': labels,
        'sm': sm,
        'bf': bf,
        'bf_unc': bf_unc,
        'cor': cor,
        'cov': cov,
        'scaling': [scaling]
    }


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

# Read in the data for each channel
channel_data = {
    "hgg": ReadHgg(),
    "wg": ReadWg(),
    "singlet": ReadSingleT()
}

# Combine the following channels
combine_channels = ["hgg", "singlet", "wg"]

# Merge the lists of labels, best-fit values, uncertainties and covariance matrices
labels = list()
for X in combine_channels:
    labels.extend(channel_data[X]['labels'])
bf = np.concatenate([channel_data[X]['bf'] for X in combine_channels])
bf_unc = np.concatenate([channel_data[X]['bf_unc'] for X in combine_channels])
cov = MergeCov([channel_data[X]['cov'] for X in combine_channels])

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

# Loop through channels in the combination, then loop through the EFT2Obs
# scaling data for each one
for X in combine_channels:
    scalings = channel_data[X]['scaling']
    for sc in scalings:
        # Loop through scaling data (most channels only have one, but e.g. Wgamma has it split in three)
        for par in sc["parameters"]:
            # Loops through the EFT coefficients defined for this parameterisation
            if par not in POIs:
                range_lo = -1.
                range_hi = 1.
                # Override with our custom ranges above
                if par in scan_ranges:
                    range_lo = scan_ranges[par][0]
                    range_hi = scan_ranges[par][1]
                POIs.append(par)
                # Insert the RooRealVar into the worksapce
                w.factory('%s[0,%g,%g]' % (par, range_lo, range_hi))
        for ib, binterms in enumerate(sc["bins"]):
            bin_label = sc['bin_labels'][ib]
            if bin_label in skip_bins:
                print('Skipping EFT parameterisation for bin %s' % bin_label)
                continue
            # Now loop through bins in the EFT2Obs data, construct a list of terms of the form:
            # 1 + sum(A_i*c_i) + sum(B_ij*c_i*c_j)
            # We will use the RooFit factory syntax for constructing a RooAddition as a sum of 
            # products. The parser is very fussy, and every term has to be a multiplication between
            # exactly two terms. Hence the "1*1" below.
            expr_parts = ['1*1']
            for term in binterms:
                val = term[0]
                # Will be a list of one or two coeffs
                vars = term[2:]
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
