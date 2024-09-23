import os
import sys
import pandas as pd
import numpy as np

# Add the parent directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hdock_utils.seq_encoder import decode_dna

def logistic(x):
    return 1 / (1 + np.exp(-x))

def map_diff_logistic(a, b, scale_factor=20):
    diff = (b - a) / scale_factor
    mapped_value = logistic(diff)
    return mapped_value

def read_hdock_scores(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        scores = []
        for line in lines[5:8]:
            tokens = line.split()
            scores.append(float(tokens[-3]))
        average_score = sum(scores) / len(scores)
        return average_score

def process_out_files(directory):
    data = []

    for filename in os.listdir(directory):
        if filename.endswith('.out'):
            aptamer_name = filename.split('-')[0]
            
            file_1c3d = f"{aptamer_name}-1c3d-Hdock.out"
            file_2a73 = f"{aptamer_name}-2a73-Hdock.out"
            
            if not (os.path.exists(os.path.join(directory, file_1c3d)) and 
                    os.path.exists(os.path.join(directory, file_2a73))):
                continue
            
            file_path_1c3d = os.path.join(directory, file_1c3d)
            target_score = read_hdock_scores(file_path_1c3d)
            
            file_path_2a73 = os.path.join(directory, file_2a73)
            non_target_score = read_hdock_scores(file_path_2a73)
            
            encoded_aptamer = aptamer_name
            aptamer_represented_by_TCGA = decode_dna(encoded_aptamer)
            
            data.append({
                'encoded_aptamer': encoded_aptamer,
                'aptamer_represented_by_TCGA': aptamer_represented_by_TCGA,
                'target_protein_name': '1c3d',
                'target_hdocklite_score': target_score,
                'non_target_name': '2a73',
                'non_target_hdocklite_score': non_target_score
            })

    return pd.DataFrame(data)

def calculate_selective_score(row):
    if pd.notna(row['target_hdocklite_score']) and pd.notna(row['non_target_hdocklite_score']):
        return map_diff_logistic(row['target_hdocklite_score'], row['non_target_hdocklite_score'])
    else:
        return None

def main():
    directory = '.'  # Current directory, change if needed
    df = process_out_files(directory)
    
    # Group by encoded_aptamer and merge target and non-target scores
    print(df.columns)
    df_merged = df.groupby('encoded_aptamer').agg({
        'aptamer_represented_by_TCGA': 'first',
        'target_protein_name': 'first',
        'target_hdocklite_score': 'first',
        'non_target_name': 'first',
        'non_target_hdocklite_score': 'first'
    }).reset_index()

    # Calculate selective score
    df_merged['selective_score'] = df_merged.apply(calculate_selective_score, axis=1)

    # Save to Excel
    df_merged.to_excel('hdock_analysis_results.xlsx', index=False)
    print("Analysis complete. Results saved to hdock_analysis_results.xlsx")

if __name__ == "__main__":
    main()