from seqfold import fold, dot_bracket

def predict_secondary_structure(dna_sequence, temperature=37):
    try:
        result = fold(dna_sequence, temp=temperature)
        structure = dot_bracket(dna_sequence, result)
        return structure
    except Exception as e:
        print(f"Error occurred: {e}")
        return None

def ssDNA_secondary_structure_prediction_mfold(dna_sequence):
    dna_sequence = dna_sequence.strip()
    temperature = 37
    result = predict_secondary_structure(dna_sequence, temperature)
    if result:
        return result
    else:
        print("无法预测二级结构。请检查输入和seqfold安装。")

# print(ssDNA_secondary_structure_prediction_mfold("AGCAGCACAGAGGTCAGATGATGGAGGTTGGTCGGGTGGGCAATCATTCTCCTATGCGTGCTACCGTGAA"))