from inference import modelInference

'''
Takes as input path, which is a path to a txt file in the format ["sequence", "label] \ n 
Reuturns two values. Seqs is a list of the sequences and labls is a list of the labels. 
They are lists and the index's correspond.
'''
def readTests(path):
    f = open(path)
    cur = f.readline() 
    seqs = []
    labls = []
    while len(cur) != 0: 
        b = cur.split('["', 2)[1]
        sequence = b.split('",')[0]
        labl = b.split(', "', 2)[1].split('"]\n')[0]
        seqs.append(sequence)
        labls.append(labl)
        cur = f.readline()
    return seqs, labls


'''
reuturns number correct and the total number of tests
'''
def runTests(seqs, labs):
    correct= 0
    for i in range (len(seqs)):
        clasif, prob = modelInference(seqs[i])
        print("progress: " + str(i/len(seqs)))
        if clasif == labs[i]:
            correct+=1
    print("achieved " + str(correct) + " correct out of " + str(len(seqs)))
    return correct, len(seqs)



if __name__ == "__main__":
    path = "E2E_test_cases.txt"
    path = "E2E_test_cases_small.txt"
    seq, labs = readTests(path)
    cor, tot = runTests(seq, labs)
    print(cor/tot)