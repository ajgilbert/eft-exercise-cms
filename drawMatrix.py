import sys
import os
import json
import ROOT
import python.plotting as plot
import numpy as np
from argparse import ArgumentParser
from python.linalg import ArrayToTMatrix


def TMatrixToTH2(matr, xlabels=[], ylabels=None, drawDiagonal=False):
    # make sure there is no memory leak when calling this function several times
    hist_name, hist_counter = matr.GetName(), 0
    while ROOT.gROOT.FindObjectAny(hist_name+str(hist_counter)): 
        hist_counter+=1
    
    Ncol, Mrow = matr.GetNcols(), matr.GetNrows()
    res = ROOT.TH2D(matr.GetName()+str(hist_counter), matr.GetTitle(), 
               Ncol, 0, Ncol, Mrow, 0, Mrow)

    # set bin labels
    if ylabels is None: 
        ylabels = xlabels[:]
    for i,xlab in enumerate(xlabels): 
        res.GetXaxis().SetBinLabel(i+1, str(xlab))
    for j,ylab in enumerate(ylabels): 
        res.GetYaxis().SetBinLabel(Mrow-j, str(ylab))
    
    # fill TH2
    for icol in range(1, Ncol+1):
        for jrow in range(1, Mrow+1-(icol-1) if drawDiagonal else Mrow+1):
            res.SetBinContent(icol, jrow, matr[Mrow-jrow][icol-1])
    
    return res


def makeMatrixPlot(matr, xlabels, ylabels, title, file_name):

    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    plot.ModTDRStyle()
    
    outdir = os.path.split(file_name)[0]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    plot.SetCorrMatrixPalette()
    #ROOT.gStyle.SetCanvasDefW(670)
    ROOT.gStyle.SetPadRightMargin(0.12)
    ROOT.gStyle.SetPaintTextFormat('.2f')
    
    canv = ROOT.TCanvas('canv', 'canv', 1000, 1000)
    hist = TMatrixToTH2(matr, xlabels, ylabels)
    hist.SetTitle(title)
    hist.SetMarkerSize(1.0)
    hist.GetXaxis().SetLabelSize(0.025)
    hist.GetYaxis().SetLabelSize(0.025)
    hist.Draw('colz text')

    latex = ROOT.TLatex()
    plot.Set(latex, NDC=None, TextFont=42, TextSize=0.03)
    plot.DrawTitle(canv, title.split(';')[0], 3)  # top right

    canv.Print(file_name)




if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--input', '-i', 
                        help='{path to JSON file}:{name of matrix}')
    parser.add_argument('--output', '-o', 
                        help='name of output file')
    parser.add_argument('--xlabels', '-x', default=[], nargs='+', 
                        help='either bin labels separated by spaces or name of bin label list in JSON file')
    parser.add_argument('--ylabels', '-y', default=[], nargs='+', 
                        help='either bin labels separated by spaces or name of bin label list in JSON file')
    parser.add_argument('--title', default='', 
                        help='title;xaxis;yaxis')
    args = parser.parse_args()

    if not args.input or ':' not in args.input:
        sys.exit('\nUsage:\n$ python drawMatrix.py '+  
                 '--input {path to JSON file}:{name of matrix}\n')

    input_file, matrix = args.input.split(':')

    if not input_file.split('.')[-1] == 'json': 
        input_file += '.json'

    if not os.path.isfile(input_file):
        sys.exit('Input file %s not found' % input_file)

    with open(input_file, 'r') as f: 
        input_dict = json.load(f)
        print('>> Reading %s' % input_file)

    if not input_dict.has_key(matrix):
        sys.exit('Matrix "%s" not found in %s' % (matrix, input_file))

    if len(args.xlabels) == 1 and input_dict.has_key(args.xlabels[0]):
        args.xlabels = input_dict[args.xlabels[0]]

    if len(args.ylabels) == 1 and input_dict.has_key(args.ylabels[0]):
        args.ylabels = input_dict[args.ylabels[0]]
    elif len(args.ylabels) == 0:
        args.ylabels = args.xlabels[:]

    if not args.output:
        output_dir, output_file = os.path.split(input_file)[0], '%s.pdf' % matrix
    else:
        output_dir, output_file = os.path.split(args.output)
        if output_dir == '':
            output_dir = os.path.split(input_file)[0]
    if not '.' in output_file:
        output_file += '.pdf'
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    makeMatrixPlot(ArrayToTMatrix(input_dict[matrix]), 
                   args.xlabels, args.ylabels, args.title, 
                   os.path.join(output_dir, output_file))



