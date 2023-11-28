import os
import numpy
from numpy import random
import scipy
import matplotlib 
import pickle
matplotlib.use('agg')
from matplotlib import pyplot as plt

import pandas as pd

import torch
import torch.nn as nn
import torchvision
import time


#Set hyper params 


#load dataset 
class genCustumDataset(Dataset):
    def __init__(self, data):
        df = pd.read_csv('data/CDH1.csv')




#construct Dataloaders?


#Evaluate model?

#Make a model 


#Train the mdoel 


