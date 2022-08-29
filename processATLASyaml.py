import yaml
import argparse
from tools import Measurement

#The argument to --measurement_file is the name of the ATLAS yaml file in the measurements directory, without the .yaml part

parser = argparse.ArgumentParser()
parser.add_argument('--measurement_file', default='')
args = parser.parse_args()


measurement = Measurement.fromYAML('measurements/{}.yaml'.format(args.measurement_file))
for i in range(measurement.nbins):
    measurement.bf[i]= measurement.bf[i]/measurement.sm[i]
    for j in range(measurement.nbins):
        if j==i:
            measurement.cov[i][j] = measurement.cov[i][j]/(measurement.sm[i]*measurement.sm[j])
        else:
            measurement.cov[i][j] = 0.

measurement.writeToYAML('measurements/{}_parsed_tmp.yaml'.format(args.measurement_file))
