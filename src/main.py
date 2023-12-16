import os
import numpy as np
import matplotlib 
matplotlib.use('agg')
from matplotlib import pyplot as plt

import pandas as pd

import torch
import torch.nn as nn
from torch.utils.data import Dataset

pd.options.mode.chained_assignment = None # Disable warnings

import sys 
np.set_printoptions(threshold=sys.maxsize) #numpy printing settings

#Set hyper params 
inSize = 20
outSize = 5
learning_rate = 1e-5
batch_size = 512
epochs = 20

##gloabl variable data
label_name_ref = []
gene_name_ref = []
dictRef = {'A':1, 'C':2, 'G':3, 'T':4}
# device = (
#     "cuda"
#     if torch.cuda.is_available()
#     else "mps"
#     if torch.backends.mps.is_available()
#     else "cpu"
# )
device = "cpu" #coding CPU for now - should use cuda on HPC, but local mps is slow


def dropIndixes(df, col):
    inds_to_drop= []
    for i, v in enumerate(df[col]):
        # print(v)
        if type(v) == float:
            inds_to_drop.append(i)
    for i, v in enumerate(df[col]):
        # print(str(i) + " length of v: " + str(v))
        if (i not in inds_to_drop):
            if (len(v) != 1 or v[0] == '-'):
                inds_to_drop.append(i)
    return inds_to_drop

def numberizeCol(df, col):
    for v in (df): 
        v[col] = dictRef[v[col]]
        # print(i)
        # print(v)

    return df


#load dataset 
class genCustDataset(Dataset):
    def __init__(self, data_path):
        df = pd.DataFrame()
        
        for dPath in data_path:
            # print(dPath)
            df_cur = pd.read_csv('data/' + dPath)
            df = pd.concat([df, df_cur], ignore_index=True)
        # df = pd.read_csv(data_path) #This should be implemented to read thru a list: lambda x: for i in x: i (data_path)
        
        #Pre processing
        inds_to_drop= []
        #Remove indices from 
        inds_to_drop.extend(dropIndixes(df, 'Ref'))
        inds_to_drop.extend(dropIndixes(df, 'Var'))

        df = df.drop(inds_to_drop)
        # print(df)
        df = df.drop(columns = ['Variant Type', 'Allele Freq (T)', 'Ref Start']) #'Reference',
        # df = df.drop(columns = ['Cancer Type Detailed', 'Variant Type', 'Reference', 'Ref Start'])
        # print(df)
        # print(df)

        #label preprocessing
        label = "Cancer Type Detailed"
        data = df.loc[:, df.columns != label]
        labels = df[label].values
        for i, v in enumerate(labels): 
            if v in label_name_ref:
                labels[i]= label_name_ref.index(v)
            else: 
                labels[i] = len(label_name_ref)
                label_name_ref.append(v)
        df = df.drop(columns = ['Cancer Type Detailed'])
        # print(df)

        #convert ACGT to numbers
        # print(df.values[0])
        mat = numberizeCol(df.values, df.columns.get_loc('Ref'))
        mat = numberizeCol(mat, df.columns.get_loc('Var'))

        #Gene name numberization

        geneCol = df.columns.get_loc('Gene Name')
        for i, v in enumerate(mat): 
            curV= v[geneCol]
            if curV in gene_name_ref:
                v[geneCol]= gene_name_ref.index(curV)
            else: 
                v[geneCol] = len(gene_name_ref)

                gene_name_ref.append(curV)

        
        #make custum dataset
        # print(mat)
        mat = mat/(mat.sum(axis=0))
        # print(labels)
        # print(type(labels))
        labels_onehot = np.zeros((labels.size, len(label_name_ref)))
        labels_onehot[np.arange(labels.size), labels.astype(int)] = 1

        # print(labels_onehot)
        self.x_train=torch.tensor(mat.astype(float), dtype=torch.float32)
        self.y_train=torch.tensor(labels_onehot.astype(float), dtype=torch.float32)
        
    def __len__(self):
        return len(self.y_train)
    # def len(self):
    #     return len(self.y_train)
    def __getitem__(self, ind):
        return (self.x_train[ind], self.y_train[ind])

#Returns dataloaders of size 0.8, 0.1, and 0.1 for train, test, and validation respectively. Sizes are percents of length of entire dataset
#This is taking the custom dataset and splitting it randomly into a training and test/valid sets 
# Also returning dataloaders, which are ideal to use internally within pytorch training
def getDataloaders(data): 
    dSets = genCustDataset(data) 
    train_size = 0.8
    test_size = 0.1
    validation_size = 0.1

    tr_dlo, te_dlo, va_dlo = torch.utils.data.random_split(dSets, [train_size, test_size, validation_size])
    # print(type(tr_dlo))

    
    tr_dload = torch.utils.data.DataLoader(tr_dlo, batch_size=batch_size, shuffle=True)
    te_dload= torch.utils.data.DataLoader(te_dlo, batch_size=100, shuffle=False)
    va_dload= torch.utils.data.DataLoader(te_dlo, batch_size=100, shuffle=False)

    return tr_dload, te_dload, va_dload



#Model 
class NeuralNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(6, 512),
            nn.BatchNorm1d(512),
            nn.LeakyReLU(0.01),
            # nn.Sigmoid(),
            nn.Linear(512, 4096),
            nn.BatchNorm1d(4096),
            nn.LeakyReLU(0.01),
            # nn.Linear(8192, 16348),
            # nn.BatchNorm1d(16348),
            # nn.LeakyReLU(0.01),
            # nn.Sigmoid(),
            nn.Linear(4096, len(label_name_ref)), #confirm the outsize is accurate on runtime 
            nn.Softmax(dim=1),
        )


    def forward(self, x):
        # x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits


def train(tr, te, model, loss_fn, optim):
    model.train()
    # loss = torch.zeros(1, requires_grad=True)
    avLoss = 0
    for i, (data, label) in enumerate(tr):
        pred = model(data)
        # print(pred)
        # print(label)
        # print(pred.argmax(1).type(torch.FloatTensor))
        # print(label.type(torch.FloatTensor))
        #pred.argmax(1).to(torch.float32)
        loss = loss_fn(pred, label)
        # print(loss)
        avLoss += loss
        optim.zero_grad()       
        loss.backward()
        optim.step()
        
        # print(i)
        # if i % 100 == 0: #UPDATE
        #     loss, current = loss.item(), (i + 1) * len(data)
        #     print(f"loss: {loss:>5f}  [{current:>2d}/{len(tr):>2d}]")
    # preds = model(tr)
    return avLoss/len(tr)

#UDPATE 
def test_loop(te, model, loss_fn):
    model.eval()
    size = 0
    num_batches = len(te)
    test_loss, correct = 0, 0

    with torch.no_grad():
        for X, y in te:
            pred = model(X)
            # print(X)
            # print(pred)
            # print(pred.argmax())
            # print(y)
            # print((pred.argmax(0) == y).type(torch.float).sum().item())
            test_loss += loss_fn(pred, y).item()
            correct += (pred.argmax(1) == y.argmax(1)).type(torch.float).sum().item()
            # print(pred.argmax() == y)
            # print(correct)
            # print(pred)

            # print(pred)
            # print(y)
            size += len(y)

    test_loss /= num_batches
    correct /= size
    print(f"Test Error: \n Accuracy: {(100*correct):>0.5f}%, Avg loss: {test_loss:>8f} \n")
    return test_loss, correct



if __name__ == '__main__':
    data_path = 'data/'
    # paths = os.listdir(data_path)
    
    
    paths = ['AllTP53.csv', 'AllGATA3.csv', 'AllCDH1.csv', 'AllPIK3CA.csv']
    # print(paths)
    # paths = ['data/CDH1.csv', 'data/CDH1.csv', 'data/CDH1.csv', 'data/CDH1.csv']

    tr, te, va = getDataloaders(paths)
    # print(tr)

    model = NeuralNetwork().to(device)
    # Initialize the loss function
    loss_fn = nn.CrossEntropyLoss()

    optim = torch.optim.SGD(model.parameters(), lr=learning_rate, momentum=0.97)



    # print(label_name_ref)
    print(tr)
    # print(len(te))
    loss_t = []
    corr_t = []
    loss_tr = []

    for iter in range (epochs):
        print("Iter: "+ str(iter))
        lossTr = train(tr, te, model, loss_fn, optim)
        loss, corr = test_loop(te, model, loss_fn)
        loss_t.append(loss)
        corr_t.append(corr)
        loss_tr.append(lossTr.item())

    print("Model wieghts")

    #save model
    path = 'CancerIdentifier.pth'
    torch.save(model.state_dict(), path)


    # plt.plot(torch.arange(epochs), loss_t, '-c', label='Approximated Testing Loss')
    # plt.title("Loss")
    # plt.legend(loc='upper right')
    # plt.xlabel("Epochs")
    # plt.ylabel("Loss")
    # # plt.savefig('Loss.png')
    # # plt.clf()
    # plt.show()

    # plt.plot(torch.arange(epochs), corr_t, '-c', label='Approximated Testing Accuracy')
    # plt.title("Accuracy")
    # plt.legend(loc='upper right')
    # plt.xlabel("Epochs")
    # plt.ylabel("Accuracy")
    # # plt.savefig('Acc.png')
    # # plt.clf()
    # plt.show()

    # plt.plot(torch.arange(epochs), loss_tr, '-c', label='Approximated Training Loss')
    # plt.title("Loss")
    # plt.legend(loc='upper right')
    # plt.xlabel("Epochs")
    # plt.ylabel("Loss")
    # # plt.savefig('LosTraoining.png')
    # # plt.clf()
    # plt.show()









    # for x, y in va:
    #     preds = model(x)


#Evaluate model

#Make a model 


#Train the mdoel 


