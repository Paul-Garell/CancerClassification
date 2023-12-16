from seq import align


# list of [Chromosome, Start Pos, End pos, Ref, Var, Gene name]
# for example [17,7578190, 7578190.0, T, C, 0.26, TP53]


# Paul inferences on model and outputs predicted cancer type probabilities
def modelInference(sequence):
    SNPs = align(sequence)
    return None




if __name__ == "__main__":
    # for testing purposes
    input = ""
    output = align(input)