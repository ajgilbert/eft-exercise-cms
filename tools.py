import numpy as np
import json
from math import sqrt

class Measurement(object):

    def __init__(self, nbins, bin_labels, sm, bf, cov):
        self.nbins = int(nbins)
        self.bin_labels = list(bin_labels)
        self.sm = np.array(sm)
        self.bf = np.array(bf)
        self.cov = np.array(cov)

    @classmethod
    def fromJSON(cls, filename):
        with open(filename) as jsonfile:
            input = json.load(jsonfile)
        return cls.fromDict(input)

    @classmethod
    def fromDict(cls, d):
        return cls(nbins=d["nbins"], bin_labels=d["bin_labels"], sm=np.array(d["sm"]), bf=np.array(d["bf"]), cov=np.array(d["cov"]))

    def writeToJSON(self, filename):
        with open(filename, 'w') as outfile:
            res = {
                "nbins": int(self.nbins),
                "bin_labels": self.bin_labels,
                "sm": self.sm.tolist(),
                "bf": self.bf.tolist(),
                "cov": self.cov.tolist(),
            }
            outfile.write(json.dumps(res, sort_keys=False, indent=2))


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
