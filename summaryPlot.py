#!/usr/bin/env python
import ROOT
# import math
import json
import argparse
import python.plotting as plot
import fnmatch
from collections import OrderedDict

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(0)

default_bar_styles = {
    '2sig_Error': {
        'LineWidth': 2,
        'LineColor': ROOT.kBlack,
        'MarkerSize': 0
    },
    '2sig_OtherLimit': {
        'LineWidth': 2,
        'LineColor': ROOT.kBlack,
        'MarkerSize': 0
    },
    'Error': {
        #'LineWidth': 2,
        'LineWidth': 4,
        'LineColor': ROOT.kBlack,
        'MarkerSize': 0
    },
    'OtherLimit': {
        'LineWidth': 4,
        'LineColor': ROOT.kBlack,
        'MarkerSize': 0
    },
    'Stat': {
        'LineWidth': 12,
        'LineColor': ROOT.kAzure+8,
        'MarkerSize': 0
    },
    'Syst': {
        'LineWidth': 20,
        'LineColor': ROOT.kRed-4,
        'MarkerSize': 0
    },
    'BestFit': {
        'MarkerSize': 1.1
    },
    'fixedOtherPOIError': {
        'LineWidth': 10,
        'LineColor': ROOT.kRed,
        'MarkerSize': 0
    },
    'fixedOtherPOI2sig_Error': {
        'LineWidth': 6,
        'LineColor': ROOT.kRed,
        'MarkerSize': 0
    }
}

default_bar_labels = {
    'fixedOtherPOI2sig_Error': '',
    'fixedOtherPOIError': '',
    'Error': '#pm1#sigma (stat #oplus syst)',
    '2sig_Error': '#pm2#sigma (stat #oplus syst)',
    'Error': '#pm1#sigma (stat #oplus syst)',
    'Stat': '#pm1#sigma (stat)',
    'Syst': '#pm1#sigma (syst)'
}


def Translate(name, ndict):
    return ndict[name] if name in ndict else name


def ParseDictArgs(str):
    return {x.split('=')[0]: eval(x.split('=')[1]) for x in str.split(',')}


def LoadTranslations(jsonfilename):
    with open(jsonfilename) as jsonfile:
        return json.load(jsonfile)


def GetListOfPOIsFromJsonFile(jsonfilename='observed.json', model='A1_mu', select=None):
    res = []
    with open(jsonfilename) as jsonfile:
        full = json.load(jsonfile, object_pairs_hook=OrderedDict)[model]
    if select is None:
        return list(full.keys())
    elif '*' in select:
        for key in reversed(full.keys()):
            if (fnmatch.fnmatch(key, select)):
                res.append(key)
    elif ',' in select:
        for key in select.split(','):
            if key in full:
                res.append(key)
    else:
        res.append(key)
    return res


def CopyDataFromJsonFile(jsonfilename='observed.json', model='A1_mu', pois=[]):
    res = []
    with open(jsonfilename) as jsonfile:
        full = json.load(jsonfile)[model]
        for poi in pois:
            res.append(dict(full[poi]))
            res[-1]['Name'] = poi
    return res


def MakeFrame(npois, xmin, xmax, xtitle='Parameter value', frac=0.8):
    hframe = ROOT.TH2F("hframe", "hframe", 6, xmin, xmax, npois+1, 0, npois+1)
    hframe.GetYaxis().SetLabelSize(0)
    hframe.GetYaxis().SetTickLength(0)
    hframe.GetXaxis().SetTitle(xtitle)
    hframe.frac = frac
    return hframe


def YPos(i, npois, hframe):
    ymin = hframe.GetYaxis().GetXmin()
    ymax = hframe.GetYaxis().GetXmax()

    height = (ymax - ymin) * hframe.frac
    per_entry = height / float(npois)
    return ymin + float(i + 0.5) * per_entry


def YEntryHeight(npois, hframe):
    ymin = hframe.GetYaxis().GetXmin()
    ymax = hframe.GetYaxis().GetXmax()
    height = (ymax - ymin) * hframe.frac
    return height / float(npois)


def MakeGraph(drawlist, hframe, label='Error', valid_checks=[]):
    gr_bar = ROOT.TGraphAsymmErrors(len(drawlist))
    for i, info in enumerate(drawlist):
        ypos = YPos(i, len(drawlist), hframe)
        if 'OtherLimit' in label:
            gr_bar.SetPoint(i, (info["%sLo" % label] + info["%sHi" % label]) / 2., ypos)
            err_lo = (gr_bar.GetX()[i] - info["%sLo" % label])
            err_hi = (info["%sHi" % label] - gr_bar.GetX()[i])
        else:
            gr_bar.SetPoint(i, info["Val"], ypos)
            err_lo = -1.0 * info["%sLo" % label]
            err_hi = info["%sHi" % label]
        valid_lo = True
        valid_hi = True
        for v in valid_checks:
            if not info["%sLo" % v]:
                valid_lo = False
            if not info["%sHi" % v]:
                valid_hi = False
        if not valid_lo:
            print 'Warning, entry %sLo for %s is not valid' % (label, info["Name"])
            # err_lo = 0.0
        if not valid_hi:
            print 'Warning, Entry %sHi for %s is not valid' % (label, info["Name"])
            # err_hi = 0.0
        gr_bar.SetPointError(i, err_lo, err_hi, 0., 0.)
    return gr_bar


def FindInvalidRegions(drawlist, xmin, xmax, valid_checks=[]):
    invalid_regions = []
    for i, info in enumerate(drawlist):
        valid_lo = True
        valid_hi = True
        for v in valid_checks:
            if not info["%sLo" % v]:
                valid_lo = False
            if not info["%sHi" % v]:
                valid_hi = False
        if not valid_lo:
            invalid_regions.append((xmin, info["Val"] + info["%sLo" % v.replace('Valid', '')], i))
        if not valid_hi:
            invalid_regions.append((info["Val"] + info["%sHi" % v.replace('Valid', '')], xmax, i))
    return invalid_regions


def MakeBestFitGraph(drawlist, hframe, label='Val'):
    gr_fit = ROOT.TGraph(len(drawlist))
    for i, info in enumerate(drawlist):
        ypos = YPos(i, len(drawlist), hframe)
        gr_fit.SetPoint(i, info[label], ypos)
    return gr_fit


def MakeYaxis(N, hframe, bin_labels=[], label_size=1.0):
    frac = hframe.frac
    ymin = hframe.GetYaxis().GetXmin()
    ymax = hframe.GetYaxis().GetXmax()
    xmin = hframe.GetXaxis().GetXmin()
    height = (ymax - ymin) * frac
    gaxis = ROOT.TGaxis(xmin, 0, xmin, ymin + height,
                   0, N, N, '-M')
    for i in xrange(N):
        gaxis.ChangeLabel(i + 1, -1, -1, -1, -1, -1, ROOT.TString(bin_labels[i]))
    gaxis.SetLabelFont(42)
    gaxis.SetLabelSize(gaxis.GetLabelSize() * label_size)
    return gaxis


def MakeLegend(pad, xlo, xhi, yhi, topfrac=0.2):
    frame_h = 1. - pad.GetBottomMargin() - pad.GetTopMargin()
    frame_frac = pad.GetBottomMargin() + (frame_h * (1. - topfrac))
    return ROOT.TLegend(xlo, frame_frac, xhi, yhi, '', 'NBNDC')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', nargs='*', help='input json file')
    parser.add_argument('--template', default='default', help='input json file')
    parser.add_argument('--output', '-o', help='name of the output file to create')
    parser.add_argument('--translate', '-t', help='JSON file for remapping of parameter names')
    parser.add_argument('--show-bars', default='Syst,Stat,2sig_Error,Error', help='JSON file for remapping of parameter names')
    parser.add_argument('--legend', default='Error,Syst,Stat,2sig_Error', help='JSON file for remapping of parameter names')
    parser.add_argument('--vlines', nargs='*', help='JSON file for remapping of parameter names')
    parser.add_argument('--hlines', nargs='*', help='JSON file for remapping of parameter names')
    parser.add_argument('--texlabels',nargs='*', help='Additional labels on the canvas')
    parser.add_argument('--require-valid', default='2sig_ValidError,ValidError', help='JSON file for remapping of parameter names')
    parser.add_argument('--cms-label', default='Internal', help='Label next to the CMS logo')
    parser.add_argument('--height', type=int, default=700, help='Canvas height in pixels')
    parser.add_argument('--width', type=int, default=600, help='Canvas width in pixels')
    parser.add_argument('--labels', default=None, help='Label next to the CMS logo')
    parser.add_argument('--x-title', default='Parameter value', help='Label next to the CMS logo')
    parser.add_argument('--x-range', default='0,10', help='Label next to the CMS logo')
    parser.add_argument('--left-margin', default=0.2, type=float, help='Left pad margin')
    parser.add_argument('--subline', default='35.9 fb^{-1} (13 TeV)', help='Label next to the CMS logo')
    parser.add_argument('--extra-text', nargs='*', help='Text:SIZE:X:Y')
    parser.add_argument('--frame-frac', type=float, default=0.8, help='Fraction of the frame y height the graphs will occupy')
    parser.add_argument('--table', default=None, help='Draw table of numeric values, with opts SIZE')
    args = parser.parse_args()

    # Dictionary to translate parameter names
    translate = {} if args.translate is None else LoadTranslations(args.translate)

    # Set the global plotting style
    if args.template == 'A1_5PD':
        plot.ModTDRStyle(l=args.left_margin, b=0.123, height=args.height, width=args.width, t=0.05) # for prod x decay plot
    else:
        plot.ModTDRStyle(l=args.left_margin, b=0.10, height=args.height, width=args.width, t=0.05)
    ROOT.gStyle.SetNdivisions(510, 'XYZ')

    canv = ROOT.TCanvas(args.output, args.output)
    pads = plot.OnePad()

    poilist = []
    drawlist = []

    # First build a list of the POIs we want to draw
    for i, iargs in enumerate(args.input):
        jsonfilename = iargs.split(':')[0]
        model, pattern = iargs.split(':')[1].split('/')
        pois = GetListOfPOIsFromJsonFile(jsonfilename, model, pattern)
        poilist.extend(pois)
        drawlist.extend(CopyDataFromJsonFile(jsonfilename, model, pois))

    print drawlist

    N = len(poilist)

    xmin = float(args.x_range.split(',')[0])
    xmax = float(args.x_range.split(',')[1])

    hframe = MakeFrame(N, xmin, xmax, xtitle=args.x_title, frac=args.frame_frac)

    relabel = {}
    if args.labels is not None:
        labels = args.labels.split(',')
        for l in labels:
            relabel[int(l.split('=')[0])] = l.split('=')[1]

    hframe.Draw()

    bin_labels = [Translate(poi, translate) for poi in poilist]
    for i, new_label in relabel.iteritems():
        bin_labels[i] = new_label

    gaxis = MakeYaxis(N, hframe, bin_labels=bin_labels, label_size=1.0)
    gaxis.Draw()

    for vline in args.vlines:
        xlines = vline.split(':')[0].split(',')
        linestyle = ParseDictArgs(vline.split(':')[1])
        line = ROOT.TLine()
        plot.Set(line, **linestyle)
        for x in xlines:
            if args.template == 'A1_5PD':
                line.DrawLine(float(x), 0., float(x), float(N) * YEntryHeight(N, hframe))
                # line.DrawLine(float(x), 0., float(x), float(N)-0.7) # for prod x decay
            else:
                line.DrawLine(float(x), 0., float(x), float(N) * YEntryHeight(N, hframe))


    graphs = []
    bars = args.show_bars.split(',')

    valid_checks = [X for X in args.require_valid.split(',') if X != '']

    for bar in bars:
        gr_bar = MakeGraph(drawlist, hframe, bar, valid_checks=valid_checks)
        plot.Set(gr_bar, **default_bar_styles[bar])
        graphs.append(gr_bar)


    gr_fit = MakeBestFitGraph(drawlist, hframe)
    plot.Set(gr_fit, **default_bar_styles['BestFit'])

    ROOT.gStyle.SetHatchesLineWidth(1)
    ROOT.gStyle.SetHatchesSpacing(3)

    invalid_regions = FindInvalidRegions(drawlist, xmin, xmax, valid_checks)

    if len(invalid_regions):
        gr_invalid = ROOT.TGraphAsymmErrors(len(invalid_regions))
        for i, inv in enumerate(invalid_regions):
            gr_invalid.SetPoint(i, inv[0], YPos(inv[2], N, hframe))
            gr_invalid.SetPointError(i, 0.0, inv[1] - inv[0], 0.5 * YEntryHeight(N, hframe), 0.5 * YEntryHeight(N, hframe))

        gr_invalid.SetFillColor(17)
        gr_invalid.SetFillStyle(3144)
        gr_invalid.Draw('2SAME')

    for gr in graphs:
        gr.Draw('ZPSAME')
    gr_fit.Draw('PSAME')

    pads[0].RedrawAxis()

    plot.DrawTitle(pads[0], args.subline, 3)

    if args.template == 'A1_5PD':
        legend = MakeLegend(pads[0], xlo=0.53, xhi=0.95, yhi=0.945) #for prod x decay
    else:
        legend = MakeLegend(pads[0], xlo=0.66, xhi=0.95, yhi=0.945)
    legend.SetFillStyle(0)

    legend.AddEntry(gr_fit, 'Observed', 'P')
    for bar in args.legend.split(','):
        i = bars.index(bar)
        legend.AddEntry(graphs[i], default_bar_labels[bar], 'L')
    legend.Draw()

    plot.DrawCMSLogo(pads[0], 'CMS',
                     args.cms_label, 11, 0.045, 0.035, 1.2, '', 1.3)

    if args.table is not None:
        table_args = args.table.split(':')
        valtxt = ROOT.TLatex()
        with_statsyst = 'Stat' in bars and 'Syst' in bars
        # pavetxt = ROOT.
        plot.Set(valtxt, TextFont=42, TextSize=float(table_args[0]), TextAlign=32)
        if with_statsyst:
            xleft = (xmax - xmin) * 0.60 + xmin
            xstart = (xmax - xmin) * 0.75 + xmin
            xstart_stat = (xmax - xmin) * 0.85 + xmin
            xstart_syst = (xmax - xmin) * 0.95 + xmin
        else:
            xleft = (xmax - xmin) * 0.80 + xmin
            xstart = (xmax - xmin) * 0.95 + xmin
        box = ROOT.TBox(xleft, 0, xmax, float(N) * YEntryHeight(N, hframe) - YEntryHeight(N, hframe) * 0.05)
        plot.Set(box, FillColor=0, LineWidth=0)
        box.Draw()
        i_err = bars.index('Error')
        for i, info in enumerate(drawlist):
            valtxt.DrawLatex(xstart, float(gr_fit.GetY()[i]), '%.2f^{#plus%.2f}_{#minus%.2f}' % (info['Val'], abs(info['ErrorHi']), abs(info['ErrorLo'])))
            if with_statsyst:
                valtxt.DrawLatex(xstart_stat, float(gr_fit.GetY()[i]), '{}^{#plus%.2f}_{#minus%.2f}' % (abs(info['StatHi']), abs(info['StatLo'])))
                valtxt.DrawLatex(xstart_syst, float(gr_fit.GetY()[i]), '{}^{#plus%.2f}_{#minus%.2f}' % (abs(info['SystHi']), abs(info['SystLo'])))
        plot.Set(valtxt, TextFont=62, TextAlign=31)
        if with_statsyst:
            valtxt.DrawLatex(xstart_stat, YEntryHeight(N, hframe) * N, 'Stat')
            valtxt.DrawLatex(xstart_syst, YEntryHeight(N, hframe) * N, 'Syst')
        line = ROOT.TLine()
        # plot.Set(line, **linestyle)
        line.DrawLine(xleft, 0., xleft, float(N) * YEntryHeight(N, hframe))
        remove_x_labels = 0
        for il in xrange(remove_x_labels):
            hframe.GetXaxis().ChangeLabel(-1 * il, -1, -1, -1, -1, -1, " ")
        # n_labels = hframe.GetXaxis().GetLabels().GetSize()
        # print n_labels


    if args.hlines == None:
        args.hlines = []
    for hline in args.hlines:
        ylines = hline.split(':')[0].split(',')
        linestyle = ParseDictArgs(hline.split(':')[1])
        line = ROOT.TLine()
        plot.Set(line, **linestyle)
        for y in ylines:
            #line.DrawLine(xmin, float(y), xmax, float(y))
            line.DrawLine(xmin, float(y)*YEntryHeight(N, hframe), xmax, float(y)*YEntryHeight(N, hframe))


    if args.texlabels == None:
        args.texlabels = []
    lbltxt = ROOT.TLatex()
    lbltxt.SetTextFont(42)
    lbltxt.SetTextAngle(90)
    lbltxt.SetTextSize(0.031*1.5)
    lbltxt.SetNDC()
    for lbl in args.texlabels:
        all_labels = lbl.split(',')
        for lbltmp in all_labels:
            x,y,lab = lbltmp.split(':')
            lbltxt.DrawLatex(float(x),float(y),lab)

    lbltxt.SetTextAngle(0)
    if args.extra_text is not None:
        for extra in args.extra_text:
            extra_split = extra.split(':')
            lbltxt.SetTextSize(float(extra_split[1]))
            lbltxt.DrawLatex(float(extra_split[2]), float(extra_split[3]), extra_split[0])


    pads[0].GetFrame().Draw()
    canv.Print('.pdf')
    #canv.Print('.png')

