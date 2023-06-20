#!/usr/bin/env python

"""
:Author: Lea Burkard
:Date: 01.11.2022
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.signal import savgol_filter
import glob, os.path
import sys, os
import argparse

import matplotlib
import seaborn as sns

from scipy import stats, signal, fft
from statsmodels.tsa.filters.filtertools import recursive_filter

from spec_pgram_mod import spec_pgram



### params -> read in from config file or command line arguments 
parser = argparse.ArgumentParser()
parser.add_argument("-r","--run", dest="run", help="Which run ID to use") #run = "DELFI_breast_cancer_corrected"
parser.add_argument("-s","--score", dest="score", help="Which score to use (MIDPOINT, WPS, COV), default:MIDPOINT", default="MIDPOINT") #score = "MIDPOINT"
parser.add_argument("-a","--amplitude", dest="amplitude", help="Which method for calculating fft amplitude (FFT, spec_pgram), default:FFT", default="FFT") #amplitude = "FFT"
parser.add_argument("-p","--path", dest="scorepath", help="path to the directory where the score is calculated") #
options = parser.parse_args()

overlay_mode = "mean"
mean_norm = True
background_norm = False
rolling = False
smoothing = True
win_len = 165
poly = 3


### functions
def calculate_flanking_regions(val: int):
    """Calculates flanking regions for point of interest.

    Args:
        val (int): should be length of value vector

    Raises:
        TypeError: Only integers are allowed

    Returns:
        [iterator]: range of values around center point (e.g. range(-1000,1000))
    """

    if not isinstance(val, int):
        raise TypeError("Only integers are allowed")

    if val % 2 == 0:
        flank = int(val / 2)
        region = range(-flank, flank)
    elif val % 2 == 1:
        flank_l = int(val / 2 - 0.5)
        flank_r = int(val / 2 + 0.5)
        region = range(-flank_l, flank_r)
    return region

taps = 1/np.arange(start=5,stop=100,step=4)
def rec_fil(x):
    xtmp = x.iloc[1:301].append(x.iloc[1:])
    ftmp = recursive_filter(xtmp,taps)[300:]
    return ftmp

def calculate_FFT(indf):
    indf = indf.apply(lambda x: rec_fil(x), axis=1)
    indf = indf.apply(lambda x: x-stats.trim_mean(x,0.1), axis=1)
    specDF = indf.apply(lambda x: spec_pgram(x,pad=0.3,taper=0.3,spans=2,plot=False,detrend=True,demean=True), axis=1)
    freq = specDF.apply(lambda x: 1/x["freq"])
    header = pd.DataFrame(freq.values.tolist(), index=freq.index).apply(lambda x: np.unique(x)).iloc[0,:].astype(int)
    datadf = specDF.apply(lambda x: x["spec"])
    datadf = pd.DataFrame(datadf.values.tolist(), columns=header)
    datadf.index = indf.index
    tmpdf = datadf.loc[:,(datadf.columns >= 120) & (datadf.columns <= 240) ]
    return tmpdf


### analysis
features = pd.DataFrame()
#path = "/data/gpfs-1/groups/ag_kircher/cfDNA-analysis/lea/masterthesis/mp_cov_pipeline/results/intermediate/"+options.run+"/table/GRCh38/target/*_"+options.score+".csv.gz"
path = options.scorepath+"/results/intermediate/"+options.run+"/table/GRCh38/target/*_"+options.score+".csv.gz"
outfile = "../features/"+options.run+"_"+options.score+"_"+options.amplitude+"_features.csv"

# prepare data frame
c=0
for file in glob.glob(path):
    c=c+1
    sample_id = file.split("--")[1].split(".")[0]
    tf_id = file.split("--")[0].split("/")[-1]
    print(c,file.split("/")[-1])
    
    if overlay_mode == "mean":
        if mean_norm:
            sample_a = pd.read_csv(file, header=None).iloc[:,1:]#.mean(numeric_only=True)
            mean_value = np.nanmean(sample_a.astype(float).values)
            sample_a = sample_a / mean_value
            sample_a = sample_a.mean(numeric_only=True)
        else:
            sample_a = pd.read_csv(file, header=None).iloc[:,1:].mean(numeric_only=True)
        if background_norm:
            sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].mean(numeric_only=True, axis=1)
            sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")
    
    elif overlay_mode == "median":
        sample_a = pd.read_csv(file, header=None).iloc[:,1:].median(numeric_only=True)
        if background_norm:
            sample_b = pd.read_csv(path_b, header=None).iloc[:,1:].median(numeric_only=True, axis=1)
            sample = pd.DataFrame(sample_a / stats.trim_mean(sample_b, 0.1), columns=["value"])
        else:
            sample = pd.DataFrame(sample_a, columns=["value"])
        sample["position"] = calculate_flanking_regions(len(sample))
        sample=sample.set_index("position")

    if smoothing:
        if rolling:
            sample = sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly)) - sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly)).rolling(1000, center=True).median() 
        else:
            sample = sample.apply(lambda x:savgol_filter(x,window_length=win_len, polyorder=poly))
    else:
        if rolling:
            sample = sample - sample.rolling(1000, center=True).median(numeric_only=True)
    

    phenotype = sample_id.split("_")[1]
    phenotype_flag = True
    if phenotype == "h":
        features.loc[sample_id,"phenotype"] = 0
    elif phenotype == "c":
        features.loc[sample_id,"phenotype"] = 1
    else:
        print("Error. No classification possible. Sample is discarded.")
        phenotype_flag = False

    if phenotype_flag:
        data = sample.reset_index(drop = False)
        data['position'] = data['position'].astype(int)
        central_cov = data[(data["position"] >= -30) & (data["position"] <= 30)]
        central_cov = central_cov.set_index('position')
        features.loc[sample_id,"central_coverage_"+tf_id] = central_cov["value"].mean(skipna=False)
        
        mean_cov = data[(data["position"] >= -1000) & (data["position"] <= 1000)]
        mean_cov = mean_cov.set_index('position')
        features.loc[sample_id,"mean_coverage_"+tf_id] = mean_cov["value"].mean(skipna=False)
        
        amplitude = data[(data["position"] >= -960) & (data["position"] <= 960)]
        amplitude = amplitude.set_index('position')

        if(options.amplitude == "FFT"):
            fft_res = np.fft.fft(amplitude["value"])
            features.loc[sample_id,"amplitude190_"+tf_id] = np.abs(fft_res[10])
            #features.loc[sample_id,"amplitude160_"+tf_id] = np.abs(fft_res[12])

            arr = np.abs(fft_res)
            max_fft = arr[8:16].max() #1920/120 = 16, 1920/240 = 8
            amp = np.where(arr == max_fft)[0][0]
            features.loc[sample_id,"nucleosome_spacing_fft_"+tf_id] = round(1920/amp)

        elif(options.amplitude == "spec_pgram"):
            indf = pd.DataFrame(amplitude.reset_index(drop=True)).T
            fft_dt = calculate_FFT(indf)
            data = pd.DataFrame(fft_dt.T)
            data.columns = ["spec"]
            data = data.reset_index(drop=False)
            data.columns = ["freq","spec"]
            outdf_part = data[(data["freq"] >= 120) & (data["freq"] <= 240)]
            max_spec = outdf_part["spec"].max()
            max_freq = outdf_part.loc[outdf_part["spec"] == max_spec,"freq"].values[0]
            amp190 = outdf_part[outdf_part["freq"] == 192]["spec"].values[0]

            features.loc[sample_id,"amplitude190_"+tf_id] = amp190
            features.loc[sample_id,"nucleosome_spacing_specprgam_"+tf_id] = max_freq
        else:
            print("Error in amplitude definition.")



features.to_csv(outfile, sep='\t')
print("Script done.")