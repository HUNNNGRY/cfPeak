#!/usr/bin/env python 

import argparse
import collections
import pandas as pd
import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
#ImportError: cannot import name 'bmat' from 'scipy.sparse.sputils' #https://github.com/aleju/imgaug/issues/770
import keras
import tensorflow as tf
from tensorflow.keras import layers
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.optimizers import Adam
from tensorflow.keras.models import load_model, save_model


def load_bigwig(path):
    bw = pyBigWig.open(path)

    return bw

def load_peak_bed(path):
    #bed6 file
    peak_bed = pd.read_csv(path, sep = "\t", header = None, index_col = None, usecols=range (6)) # only import first 6 cols
    return peak_bed

def load_tx_gn(path):
    tx_gn = pd.read_csv(path, sep = "\t", header = 0, index_col = None)

    return tx_gn

def process_peak_bed(peak_bed, tx_info, bw, input_size):
    #merge tx_info
    peak_bed = pd.merge(peak_bed, tx_info, how = "left", left_on = 0, right_on = "transcript_id")
    #filter peaks
    peak_bed["mid_point"] = peak_bed.apply(lambda x : int((x[1] + x[2])/2), axis = 1)
    peak_bed_filtered = peak_bed[ (peak_bed["mid_point"] >= input_size) & (peak_bed["tx.length"]>= peak_bed["mid_point"] + input_size)]
    #get peak coverage
    peak_bed_final = pd.DataFrame(peak_bed_filtered[[0,1,2,"transcript_type","tx.length","mid_point"]].values,index=peak_bed_filtered[3].values, columns=["tx_id", "peak_start", "peak_end", "type", "tx_length", "mid_point"])
    
    def get_coverage(row):
        try:
            tx_length = bw.chroms(row["tx_id"])
            coverage = bw.values(row["tx_id"], row["peak_start"], row["peak_end"])
        except:
            print("Error: transcript ID {0} not in the bigwig file!".format(row["tx_id"]))
            return np.nan

        try:
            coverage = bw.values(row["tx_id"], row["peak_start"], row["peak_end"])
        except:
            print("Error: interval {0} to {1} not in transcript {2}!".format(row["peak_start"], row["peak_end"], row["tx_id"]))
            return np.nan

        return bw.values(row["tx_id"], row["peak_start"], row["peak_end"])

    peak_bed_final["peak_coverage"] = peak_bed_final.apply(get_coverage, axis=1)
    # Filter Rows with NaN values in peak_coverage
    peak_bed_final = peak_bed_final.dropna(subset=["peak_coverage"])
    
    peak_bed_final["input_coverage"] = peak_bed_final.apply( lambda x : list(np.nan_to_num(bw.values(x["tx_id"], x["mid_point"] - input_size, x["mid_point"] + input_size), 0)), axis=1)

    scaler = MinMaxScaler()
    peak_bed_final["scaled_coverage"] = peak_bed_final.apply(lambda x : scaler.fit_transform(np.array(x["input_coverage"]).reshape(-1,1)).reshape(1,-1)[0], axis = 1)

    return peak_bed_final

def df_to_1D_input_data(df):
    #dataframe to CNN input or predict data
    X = np.array(df["scaled_coverage"].tolist())
    y = None

    time_steps = len(X[0])
    n_features = 1
    input_shape = (time_steps, n_features)
    X = X.reshape(X.shape[0], time_steps, n_features)
    if "label" in df.columns:
        y = np.array(df['label'].tolist())
        y = keras.utils.np_utils.to_categorical(y, 2)

    return X, y

def model_predict(peak_bed_df, model): # model = None
    x_predict,_ = df_to_1D_input_data(peak_bed_df)

    # if model == None:
    model = load_model(model)
    res = model.predict(x_predict)

    return res
    
def plot_scaled_peak(peak_df, save_path=None, index_number=True):
    #plot peaks from peak dataframe
    fig_size = int(len(peak_df) ** 0.5) + 1
    
    plt.figure(figsize=(fig_size,fig_size))
    
    if index_number:
        for i in range(len(peak_df)):
            plt.subplot(fig_size,fig_size,i+1,)
            plt.plot(peak_df["input_coverage"].values[i])
            plt.text(1,1,i,color = 'r',fontsize=8)
            plt.axvline(50,color = "r",ls = "--")
            plt.xticks([])
            plt.yticks([])
    else:
        for i in range(len(peak_df)):
            plt.subplot(fig_size,fig_size,i+1,)
            plt.plot(peak_df["input_coverage"].values[i])
            plt.axvline(50,color = "r",ls = "--")
            plt.xticks([])
            plt.yticks([])
    
    if save_path != None:
        plt.savefig(save_path)



def main():
    parser = argparse.ArgumentParser(description = 'Anomaly detection of candidate peaks by CNN model')
    parser.add_argument('--bed6', '-b', required = True, help = "Input candidate peak bed6 file")
    parser.add_argument('--bigwig', '-bw', required = True, help = "Input sample bigwig file")
    parser.add_argument('--model', '-m', default = "./model/cnn_model.h5", help = "Trained CNN model file")
    parser.add_argument('--threshold', '-t', default = 0.5, help = "threshold of CNN model prediction probability, default = 0.5")
    parser.add_argument('--half_bin_size', default = 50, help = "half of bin size: frame length of image for CNN training, default = 50")
    parser.add_argument('--tx_tab', default = "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt", 
                        help = "path to tx metadat table information, should included a header named transcript_id,transcript_type,tx.length, denotes tx_id,tx_type,tx_length, respectively. default = /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt")
    parser.add_argument('--plot_type', '-p', default = "No", choices = ["No","All","Anomaly"], help = "plot the peaks(All is 900 pic at most), default = No")
    parser.add_argument('--plot_path', '-pp', default = "./peak.pdf", help = "path to save plots, default = ./peak.pdf")
    parser.add_argument('--output', '-o', required = True, help = "output bed6 to save anomalous peaks")
    args = parser.parse_args()

    input_size = int(args.half_bin_size)

    tx_gn_path = args.tx_tab # "/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt"
    tx_gn = load_tx_gn(tx_gn_path)

    peak_bed6 = load_peak_bed(args.bed6)
    bw = load_bigwig(args.bigwig)
    peak_bed = process_peak_bed(peak_bed6, tx_gn, bw, input_size)
    res = model_predict(peak_bed, model=args.model)
    anomaly_mask = res[:,1] > float(args.threshold)
    # peak_bed6[,4] = res[:,1]  # todo: change score (5th col) to prob.
    
    #change score to prob by wtw
    prob_dict = dict(zip(peak_bed.index.values, res[:,1]))
    peak_bed6[4] = peak_bed6.apply(lambda x : prob_dict[x[3]] if (x[3] in peak_bed.index.values) else -1, axis = 1) # np.nan

    peak_bed6.loc[peak_bed6[3].isin(peak_bed.loc[anomaly_mask].index.values)].to_csv(str(args.output+".removed"),header=False,index=False,sep="\t")
    peak_bed6.loc[-(peak_bed6[3].isin(peak_bed.loc[anomaly_mask].index.values))].to_csv(args.output,header=False,index=False,sep="\t")
    print("removed ratio:"+str(sum(peak_bed6[3].isin(peak_bed.loc[anomaly_mask].index.values)))+"/"+str(peak_bed6.shape[0])+"="+str(100*sum(peak_bed6[3].isin(peak_bed.loc[anomaly_mask].index.values))/peak_bed6.shape[0])+"%")
    
    if args.plot_type == "All":
        plot_scaled_peak(peak_bed.iloc[:900], args.plot_path, True)
    elif args.plot_type == "Anomaly":
        plot_scaled_peak(peak_bed.loc[anomaly_mask], args.plot_path, True)

if __name__ == "__main__":
    main()
