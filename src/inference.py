from seq import align
import main
import torch
import torch.nn as nn

# list of [Chromosome, Start Pos, End pos, Ref, Var, Gene name]
# for example [17,7578190, 7578190.0, T, C, 0.26, TP53]
type_labels= ['Breast Invasive Ductal Carcinoma', 'Breast Invasive Lobular Carcinoma', \
            'Invasive Breast Carcinoma', 'Breast Invasive Cancer, NOS', \
            'Breast Invasive Carcinoma, NOS', 'Breast Invasive Carcinosarcoma, NOS', \
            'Metaplastic Breast Cancer', 'Breast Invasive Mixed Mucinous Carcinoma',\
            'Breast', 'Paget Disease of the Nipple', 'Solid Papillary Carcinoma of the Breast']
dictRef = {'A':1, 'C':2, 'G':3, 'T':4}
gene_name_refs={'TP53':0, 'GATA3':1, 'CDH1':2, 'PIK3CA':3}

class NeuralNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(6, 512),
            nn.BatchNorm1d(512),
            nn.LeakyReLU(0.01),
            # nn.Sigmoid(),
            nn.Linear(512, 2048),
            nn.BatchNorm1d(2048),
            nn.LeakyReLU(0.1),
            nn.Linear(2048, 2048),
            nn.BatchNorm1d(2048),
            nn.LeakyReLU(0.1),
            # nn.Linear(8192, 8192),
            # nn.BatchNorm1d(8192),
            # nn.LeakyReLU(0.01),
            # nn.Sigmoid(),
            nn.Linear(2048, 11), #confirm the outsize is accurate on runtime 
            nn.Softmax(dim=1),
        )


    def forward(self, x):
        # x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits

'''
pretty printer
'''
def prettyPrint(pred, prob):
    print("Prediction: " + pred + f".  With probability: {(100*prob):>0.2f}% \n")



def preProcess(d):
    # [[10, 8115729, 8115729, 'G', 'T', 'GATA3']]
    d[3] = dictRef[d[3]]
    d[4] = dictRef[d[4]]
    d[5] = gene_name_refs[d[5]]
    return torch.tensor([d], dtype=torch.float32)


# Paul inferences on model and outputs predicted cancer type probabilities
'''
Returns the prediction and its associated probability
'''
def modelInference(sequence):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    PATH = 'CancerIdentifier.pth'
    model = main.NeuralNetwork().to(device)
    model.load_state_dict(torch.load(PATH))
    model.eval()
    test  = model(torch.tensor([[0,0,0,0,0,0]], dtype=torch.float32))
    acc = torch.zeros_like(test)
    SNPs = align(sequence)
    print(SNPs)
    #SNPs is a list of mutations
    for sn in SNPs: 
        psn = preProcess(sn)
        acc+= model(psn)
    predIndex = torch.argmax(acc)
    acc = acc/acc.sum()
    pred = type_labels[predIndex]
    probability = torch.max(acc)
    # prettyPrint(pred, probability)
    return pred, probability





if __name__ == "__main__":
    # for testing purposes
    input = ""
    sequence = align(input)

    pred, prob = modelInference(sequence)
    prettyPrint(pred, prob)