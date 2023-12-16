from seq import align
import main
import torch
import numpy as np

# list of [Chromosome, Start Pos, End pos, Ref, Var, Gene name]
# for example [17,7578190, 7578190.0, T, C, 0.26, TP53]
type_labels= ['Breast Invasive Ductal Carcinoma', 'Breast Invasive Lobular Carcinoma', \
            'Invasive Breast Carcinoma', 'Breast Invasive Cancer, NOS', \
            'Breast Invasive Carcinoma, NOS', 'Breast Invasive Carcinosarcoma, NOS', \
            'Metaplastic Breast Cancer', 'Breast Invasive Mixed Mucinous Carcinoma',\
            'Breast', 'Paget Disease of the Nipple', 'Solid Papillary Carcinoma of the Breast']

# Paul inferences on model and outputs predicted cancer type probabilities
'''
Returns the prediction and its associated probability
'''
def modelInference(sequence):
    SNPs = align(sequence)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    PATH = 'model/CancerIdentifier.pth'
    model = main.NeuralNetwork().to(device)
    model.load_state_dict(torch.load(PATH))
    model.eval()
    test  = model([0,0,0,0,0,0])
    acc = np.zeros_like(test)
    #SNPs is a list of mutations
    for sn in SNPs: 
        acc+= model(sn)
    predIndex = np.argmax(acc)
    acc = acc/acc.sum()
    pred = type_labels[predIndex]
    probability = np.max(acc)
    return pred, probability


'''
pretty printer
'''
def prettyPrint(pred, prob):
    print("Prediction: " + pred + f".  With probability: {(100*prob):>0.2f}% \n")


if __name__ == "__main__":
    # for testing purposes
    input = ""
    sequence = align(input)

    pred, prob = modelInference(sequence)
    prettyPrint(pred, prob)