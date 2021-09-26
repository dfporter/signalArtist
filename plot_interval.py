import pandas, os, sys, matplotlib, HTSeq, glob, pyBigWig, scipy, importlib, random, gffutils, re
#import pysam, collections
import seaborn as sns
import numpy as np

from typing import List, Mapping, Union
import matplotlib.pyplot as plt

# For smoothing.
from scipy.ndimage.filters import uniform_filter1d

from . import geneDrawingUtils
importlib.reload(geneDrawingUtils)

def rc(seq):
    _dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(_dict[x.upper()] for x in seq[::-1])


def motif_density_in_iv(iv, genomic_fasta, motif='GCATG', window=50):
    
    seq = genomic_fasta[iv[0]][int(iv[1]):int(iv[2])]
    if len(iv) == 4 and iv[3] == '-':
        seq = rc(seq)
        
    arr = np.zeros(len(seq))

    for m in re.finditer(motif, seq):
        arr[m.start():m.end()] = 1

    return uniform_filter1d(arr, window)

def get_coverage_in_an_interval_for_bam_file(bam_fname, interval):
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    # Get np.array of coverage depth across interval.
    b = np.sum(samfile.count_coverage(*interval[:-1]), axis=0)
    #if norm_to_per_million:
    #    b = b * 1000000/rbfox_read_numbers[os.path.basename(bam_fname)]
    samfile.close()
    return b


def get_coverage_in_an_interval_for_bigwig_file(bw_fname, interval):
    bw = pyBigWig.open(bw_fname)
    return np.nan_to_num(bw.values(*interval[:-1]))


def coverage(fname, interval, norm_to_per_million=True):
    
    if interval[3] == '-':
        fname = re.sub('_\+\.', '_-.', fname)
    else:
        fname = re.sub('_\-\.', '_+.', fname)

    if fname.split('.')[-1] == 'bam':
        return get_coverage_in_an_interval_for_bam_file(fname, interval, norm_to_per_million=norm_to_per_million)
    return get_coverage_in_an_interval_for_bigwig_file(fname, interval)


def plot_interval(
    iv, bw_fnames={}, genomic_fasta=None, window=50, db=None, RNAs=None,
    motif='GCAGCA', plot_replicates=True):

    signalArrs = {}
    for name, bw_fname_set in bw_fnames.items():
        signalArrs[name] = [np.abs(coverage(fname, iv)) for fname in bw_fname_set]
    
    averages, smoothed, smoothed_ave = {}, {}, {}
    for name, list_of_arrays in signalArrs.items():
        if len(list_of_arrays) > 0:
            averages[name] = np.sum(list_of_arrays, axis=0)/len(list_of_arrays)
            smoothed[name] = [uniform_filter1d(arr, size=window) for arr in list_of_arrays]
            smoothed_ave[name] = uniform_filter1d(averages[name], size=window)
    names = list(smoothed_ave.keys())
    
    if genomic_fasta:
        motif_density = motif_density_in_iv(iv, genomic_fasta, motif=motif, window=window)

    # Plot.
    if plot_replicates:
        fig, ax = plt.subplots(4, 1, figsize=(9, 6), gridspec_kw={'height_ratios':[2,2,0.7,0.3]})
        motif_ax = 2
        gene_ax = 3
    else:
        fig, ax = plt.subplots(3, 1, figsize=(5, 4), gridspec_kw={'height_ratios':[2,0.7,0.3]})
        motif_ax = 1
        gene_ax = 2

    # Plot averages.
    colors = 'krbgv'
    for name, aves_arr in smoothed_ave.items():
        ax[0].plot(aves_arr, c=colors[names.index(name)], alpha=0.7, lw=0.5, label=name)

    if genomic_fasta:
        ax[motif_ax].plot(motif_density, c='g', alpha=0.7, lw=0.5)
    
    # Plot replicates in the axis below.
    if plot_replicates:
        for name, list_of_arrays in smoothed.items():
            for arr in list_of_arrays:
                ax[1].plot(arr, c=colors[names.index(name)], alpha=0.7, linestyle=':', lw=0.5)
    
    # Plot any gene that might be there.
    if db is not None:
        _gene = geneDrawingUtils.plot_gene(iv, ax[gene_ax], db, RNAs)

    for _ax in ax:
        _ax.set_xticks([], [])
        _ax.set_ylabel('RPM')
        
    ax[motif_ax].set_ylabel('Motif')
    print('-----')
    ax[0].annotate(f"{iv[0]}:{iv[1]:,}-{iv[2]:,}", xy=(0.05, 0.80), xycoords='axes fraction')
    #fig.set_figwidth(12); fig.set_figheight(6)
    
    def set_ticks(_ax):
        yticklabel = _ax.get_yticks()[-2]
        _ax.set_yticks([0, _ax.get_yticks()[-2]])
        yticklabel = int(float(yticklabel) * 1000)/1000
        _ax.set_yticklabels([0, f"{yticklabel}"])
    
    set_ticks(ax[0])
    if plot_replicates:
        set_ticks(ax[1])    
    
    os.makedirs('figs/signals/', exist_ok=True)
    fig.savefig('figs/signals/' + ' '.join([str(x) for x in iv]) + f'.{_gene}.pdf')
    plt.show()
    plt.clf(); plt.close()

    
def plot_intervals(ivs, bw_fnames={}, genomic_fasta=None, window=50, n_to_plot=None, ):
    if n_to_plot is None or (n_to_plot > len(ivs)):
        n_to_plot = len(ivs)

    os.makedirs('figs/signals/', exist_ok=True)
    
    [plot_interval(ivs[n], bw_fnames=bw_fnames, genomic_fasta=genomic_fasta, window=window) for n in range(0, n_to_plot)]
    

