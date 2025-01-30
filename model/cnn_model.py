"""
This is the model construction and training demo for the simple CNN model in the exPeak pipeline.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

import keras
import tensorflow as tf
from tensorflow.keras import layers
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.optimizers import Adam

# Data preparation

# Functions
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

def plot_scaled_peak(peak_df, save_path=None, index_number=False):
    #plot peaks from peak dataframe
    fig_size = int(len(peak_df) ** 0.5) + 1
    
    plt.figure(figsize=(fig_size,fig_size))
    
    if index_number:
        for i in range(len(peak_df)):
            plt.subplot(fig_size,fig_size,i+1,)
            plt.plot(eval(peak_df["input_coverage"].values[i]))
            plt.text(1,1,i,color = 'r',fontsize=8)
            plt.axvline(50,color = "r",ls = "--")
            plt.xticks([])
            plt.yticks([])
    else:
        for i in range(len(peak_df)):
            plt.subplot(fig_size,fig_size,i+1,)
            plt.plot(eval(peak_df["input_coverage"].values[i]))
            plt.axvline(50,color = "r",ls = "--")
            plt.xticks([])
            plt.yticks([])
    
    if save_path != None:
        plt.savefig(save_path)

# Load the training data in dataframe format
# The training data is a file with the following columns:
# index(peak_x) tx_id(ENSTxxxx)   peak_start      peak_end        type(xRNA)    tx_length       mid_point       peak_coverage(list)   label   input_coverage(list)
"""The true and false peaks should be loaded depending on your data"""
T_peak["label"] = 1
F_peak["label"] = 0

scaler = MinMaxScaler()

x_df = pd.concat([T_peak, F_peak], axis = 0)
x_df["scaled_coverage"] = x_df.apply(lambda x : scaler.fit_transform(np.array(eval(x["input_coverage"])).reshape(-1,1)).reshape(1,-1)[0], axis = 1)

#stratified training and validation division
valid_ratio = 0.2
T_df = x_df[x_df["label"]==1]
F_df = x_df[x_df["label"]==0]
F_valid_number = int(valid_ratio * len(F_df))
T_valid_number = int(valid_ratio * len(T_df))

x_valid_df = pd.concat([F_df.iloc[0:F_valid_number], T_df.iloc[0:T_valid_number]], axis = 0)
x_train_df = pd.concat([F_df.iloc[F_valid_number:], T_df.iloc[T_valid_number:]], axis = 0)

x_train, y_train = df_to_1D_input_data(x_train_df)
x_valid, y_valid = df_to_1D_input_data(x_valid_df)

# Model implementation
# The parameters of the model can be adjusted according to the data
# Here is a demo model with 2 convolutional layers and 2 max pooling layers
input_shape = (100, 1)

filters = 16

model_cnn = Sequential()
model_cnn.add(Conv1D(filters, kernel_size = 3, strides = 1, padding = 'same', activation = 'relu', input_shape = input_shape))

model_cnn.add(Conv1D(filters, kernel_size = 3, strides = 1, padding = 'same', activation = 'relu', input_shape = input_shape))

model_cnn.add(MaxPooling1D(pool_size = 2, strides = 1))

model_cnn.add(Conv1D(filters*2, 3, activation = 'relu', padding = 'same'))

model_cnn.add(Conv1D(filters*4, 3, activation = 'relu', padding = 'same'))

model_cnn.add(MaxPooling1D(pool_size = 2, strides = 2))


#model_cnn.add(Dropout(0.25))
model_cnn.add(Flatten())
model_cnn.add(Dense(128,activation = 'relu'))
#model_cnn.add(Dropout(0.5))
model_cnn.add(Dense(2, activation = 'softmax'))
model_cnn.summary()

model_cnn.compile(loss = 'categorical_crossentropy', optimizer = Adam(lr = 0.0001), metrics = ['accuracy'])
hist2 = model_cnn.fit(x_train, y_train, epochs = 30, batch_size = 32,verbose = 1, validation_data = (x_valid, y_valid))

#model saving
from tensorflow.keras.models import load_model, save_model
save_model(model_cnn, './cnn_model.h5')

# Plotting the loss
y_vloss = hist2.history['val_loss']
y_loss = hist2.history['loss']
x_len = np.arange(len(y_loss))

plt.plot(x_len, y_vloss, marker='.', c='red', label="Validation-set Loss")
plt.plot(x_len, y_loss, marker='.', c='blue', label="Train-set Loss")
plt.legend()
plt.title('cnn_epoch50_loss')
#plt.savefig('loss_allpeak.svg')
