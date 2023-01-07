import yaml
import json
import numpy as np
from tools import Measurement, ReadIndependent, ReadDependent, ReadYodaString


def ReadHgg(resources):
    # Read the HIG-19-015 'minimal merging' STXS results and correlation matrix
    with open(resources['measurement'], 'r') as f:
        vals = yaml.safe_load(f)
    with open(resources['covariance'], 'r') as f:
        corr = yaml.safe_load(f)

    # Read the labels from the hepData table. However, for now we'll override these manually
    # with labels that correspond to the bin naming in the STXS Rivet ouput from EFT2Obs.
    labels = ReadIndependent(vals, col=0)
    # Some labels combine multiple STXS bins, that are merged in the analysis. These are given
    # as comments below.
    new_labels = [
        'GG2H_0J_PTH_0_10',
        'GG2H_0J_PTH_GT10',
        'GG2H_1J_PTH_0_60',
        'GG2H_1J_PTH_60_120',
        'GG2H_1J_PTH_120_200',
        'GG2H_GE2J_MJJ_0_350_PTH_0_60',
        'GG2H_GE2J_MJJ_0_350_PTH_60_120',
        'GG2H_GE2J_MJJ_0_350_PTH_120_200',
        'GG2H_PTH_200_300',
        'GG2H_PTH_300_450',
        'GG2H_PTH_GT450', # 'GG2H_PTH_450_650', 'GG2H_PTH_GT650',
        'QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25', # + 'GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25',
        'QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25', # + 'GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25',
        'QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25', # + 'GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25',
        'QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25', # + 'GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25',
        'QQ2HQQ_GE2J_MJJ_60_120',
        'QQ2HQQ_GE2J_MJJ_GT350_PTH_GT200',
        'QQ2HLNU_PTV_0_75',
        'QQ2HLNU_PTV_75_150',
        'QQ2HLNU_PTV_GT150', #'QQ2HLNU_PTV_150_250_0J',    'QQ2HLNU_PTV_150_250_GE1J', 'QQ2HLNU_PTV_GT250', 
        'QQ2HLL', # 'QQ2HLL_PTV_0_75', 'QQ2HLL_PTV_75_150', 'QQ2HLL_PTV_150_250_0J', 'QQ2HLL_PTV_150_250_GE1J', 'QQ2HLL_PTV_GT250', 'GG2HLL_FWDH', 'GG2HLL_PTV_0_75', 'GG2HLL_PTV_75_150', 'GG2HLL_PTV_150_250_0J', 'GG2HLL_PTV_150_250_GE1J', 'GG2HLL_PTV_GT250',
        'TTH_PTH_0_60',
        'TTH_PTH_60_120',
        'TTH_PTH_120_200',
        'TTH_PTH_200_300',
        'TTH_PTH_GT300',
        'TH'
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

    # cor = ROOT.TMatrixDSym(N)
    cov = np.zeros((N, N))

    # This part is analysis dependent - have to know the (i,j) index ordering
    # in the 1D array of values
    for i in range(N):
        for j in range(N):
            # cor[i][j] = corr_vals[i * N + j]
            # compute the covariance
            cov[i][j] = corr_vals[i * N + j] * bf_unc[i] * bf_unc[j]
    
    # Load the EFT2Obs scaling json
    return Measurement(nbins=N, bin_labels=labels, sm=sm, bf=bf, cov=cov)


def ReadWg(resources):
    # Read the SMP-20-005 pTgamma x |phi_f| cross section

    # The unrolling of the 2D is phi first, then pT. To match up with the cov.
    # matrix it will be easiest to covert to pT,phi unrolled instead.
    def Reshape(arr):
        # BAD: we hardcode '3' as the number of phi bins
        return arr.reshape((3, int(len(arr) / 3))).T.reshape((-1))

    with open(resources['measurement'], 'r') as f:
        vals = yaml.safe_load(f)
    with open(resources['covariance'], 'r') as f:
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

    cov = np.zeros((N, N))

    # Ordering of terms in the 1D correlation array is non-intuitive,
    # compared to Hgg  above, hence the complicated index logic:
    idx = 0
    for i in range(N):
        for j in range(N - i):
            # cor[i][N - j - 1] = corr_vals[idx]
            # cor[N - j - 1][i] = cor[i][N - j - 1] 
            cov[i][N - j - 1] = corr_vals[idx] * bf_unc[i] * bf_unc[N - j - 1]
            cov[N - j - 1][i] = cov[i][N - j - 1]
            idx += 1

    return Measurement(nbins=N, bin_labels=labels, sm=sm, bf=bf, cov=cov)


def ReadSingleT(resources):
    # Read the TOP-17-023 top pT corss section
    with open(resources['measurement'], 'r') as f:
        vals = yaml.safe_load(f)
    with open(resources['covariance'], 'r') as f:
        cov_file = yaml.safe_load(f)
    # The SM predictions are not in HEPData, have to read them 
    # from a YODA file that can not be opened with yoda.read()
    with open(resources['prediction'], 'r') as f:
        pred_str = f.read().replace('\t',' ')

    substrings = list()
    while pred_str.find('# END') != -1:
        lenline = pred_str[pred_str.find('# END'):].find('\n')+2
        substrings.append(pred_str[:pred_str.find('# END')+lenline])
        pred_str = pred_str[pred_str.find('# END')+lenline:]
    
    pred = dict([ReadYodaString(substr) for substr in substrings])
    sm = pred['Powheg4FS']['vals']

    labels = ReadIndependent(vals, col=0)
    labels = ['pt_t_bin_%i' % X for X in range(len(labels))]
    N = len(labels)
    
    bf = ReadDependent(vals, col=0)
    bf_unc = ReadDependent(vals, col=0, error=[-1], sym_errors=True)
    
    assert(N == len(sm) == len(bf) == len(bf_unc))
    bf = bf / sm
    bf_unc = bf_unc / sm

    cov_vals = ReadDependent(cov_file, col=0)
    assert((N * N) == len(cov_vals))

    cov = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            # Here we are given the covariance of the differential cross section
            # directly - so we need to convert to ratio wrt. SM
            cov[i][j] = cov_vals[i * N + j] / (sm[i] * sm[j])
        
    return Measurement(nbins=N, bin_labels=labels, sm=sm, bf=bf, cov=cov)


with open('hepdata_inputs/resources.json','r') as f:
    resources_dict = json.load(f)

# Read in the data for each channel
channel_data = {
    'CMS_hgg_STXS': ReadHgg(resources_dict['hgg']),
    'CMS_wgamma': ReadWg(resources_dict['wg']),
    'CMS_singlet': ReadSingleT(resources_dict['single_t'])
}

for label, data in channel_data.items():
    data.writeToJSON('measurements/{}.json'.format(label))
    #data.writeToYAML('measurements/{}.yaml'.format(label))
