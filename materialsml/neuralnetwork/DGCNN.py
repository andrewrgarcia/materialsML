import pandas as pd
import numpy as np

import stellargraph as sg
from stellargraph.mapper import PaddedGraphGenerator
from stellargraph.layer import DeepGraphCNN
from stellargraph import StellarGraph

# from stellargraph import datasets

from sklearn import model_selection
from IPython.display import display, HTML

from tensorflow.keras import Model
# from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense, Conv1D, MaxPool1D, Dropout, Flatten
from tensorflow.keras.losses import binary_crossentropy
import tensorflow as tf

class GraphNet():
    def __init__(self,GRAPHS):
        self.GRAPHS = GRAPHS

    def build(self):

        generator = PaddedGraphGenerator(graphs=self.GRAPHS)

        k = 35  # the number of rows for the output tensor
        layer_sizes = [32, 32, 32, 1]

        dgcnn_model = DeepGraphCNN(
            layer_sizes=layer_sizes,
            activations=["tanh", "tanh", "tanh", "tanh"],
            k=k,
            bias=False,
            generator=generator,
        )

        x_inp, x_out = dgcnn_model.in_out_tensors()



        x_out = Conv1D(filters=16, kernel_size=sum(layer_sizes), strides=sum(layer_sizes))(x_out)
        x_out = MaxPool1D(pool_size=2)(x_out)

        x_out = Conv1D(filters=32, kernel_size=5, strides=1)(x_out)

        x_out = Flatten()(x_out)

        x_out = Dense(units=128, activation="relu")(x_out)
        x_out = Dropout(rate=0.5)(x_out)

        # predictions = Dense(units=1, activation="sigmoid")(x_out)
        predictions = Dense(units=1)(x_out)

        return Model(inputs=x_inp, outputs=predictions)


