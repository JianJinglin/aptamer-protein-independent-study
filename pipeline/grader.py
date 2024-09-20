from rna_downloader import *
from hdock_utils import *
import numpy as np
class Grader:
    def __init__(self):
        self.name = ""
        self.score = 0

    def get_dummy_score(self, apta):
        random_score = np.random.uniform(-1, 1)
        return random_score

    def get_score(self, apta):
        # file_path = get_rna_composition_and_download(apta)
        file_path = get_dna_composition_and_download(apta)
        print(file_path)

        # Multi-thread
        target_score = 0
        non_target_score = 0
        with ThreadPoolExecutor() as executor:
          futures = [executor.submit(get_scores_with_hdock, file_path, './protein_pdb_files/C3dg.pdb'),
                    executor.submit(get_scores_with_hdock, file_path, './protein_pdb_files/2a73.pdb')]
          results = [future.result() for future in futures]
          target_score = results[0]
          non_target_score = results[1]

        # self.score = -(target_score - non_target_score)
        self.score = self.map_diff_logistic(target_score, non_target_score)
        # probably we can maintain a dataframe here: aptamer, target, target_score, non-target, non_target_score

        return self.score

    def logistic(self, x):
        return 1 / (1 + np.exp(-x))

    def map_diff_logistic(self, a, b, scale_factor=20):
        diff = (b - a) / scale_factor
        mapped_value = self.logistic(diff)
        return mapped_value

    def run_hdock_parallel(self, apta):
        dna_file_path = get_dna_composition_and_download(apta)
        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(get_scores_with_hdock, dna_file_path, './protein_pdb_files/C3dg.pdb'),
                        executor.submit(get_scores_with_hdock, dna_file_path, './protein_pdb_files/2a73.pdb')]
            results = [future.result() for future in futures]
            target_score = results[0]
            non_target_score = results[1]
        return results
    
    def test_one_sequence(self, aptamer_sequence):
        print(f"Testing aptamer sequence: {aptamer_sequence}")
        results = self.run_hdock_parallel(aptamer_sequence)
        print(f"Results: {results}")
        return results

# example
grader = Grader()
aptamer_sequence = "TCACCGCCGTTTTAA"  # your aptamer sequence
grader.test_one_sequence(aptamer_sequence)