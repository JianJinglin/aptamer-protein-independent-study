import os
import subprocess

def read_hdock_scores(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        scores = []
        for line in lines[5:8]:
            tokens = line.split()
            scores.append(float(tokens[-3]))
        average_score = sum(scores) / len(scores)
        return average_score

def get_scores_with_hdock(apta_path, protein_path):
    # print(f"before apta_path:{apta_path}")
    # print(f"before protein_path:{protein_path}")
    apta_name = apta_path.split('/')[-1].split('.')[0] # encoded_name
    # print(f"apta_name:{apta_name}")
    protein_name = protein_path.split('/')[-1].split('.')[0]
    # print(f"protein_name:{protein_name}")
    hdock_out_path = f"{apta_name}-{protein_name}.out"

    # apta_path = os.path.join('Aptamer_pdb', f'{apta_name}.pdb')
    # protein_path = os.path.join('protein_pdb_files', f'{protein_name}.pdb')
    new_apta_path = os.path.join('Aptamer_pdb', f'{apta_name}.pdb').replace('\\', '/')
    new_protein_path = os.path.join('protein_pdb_files', f'{protein_name}.pdb').replace('\\', '/')
    
    # print(f"new_apta_path:{new_apta_path}")
    # print(f"new_protein_path:{new_protein_path}")
    print(f"Starting hdock: {hdock_out_path}")

    # Early return (Cache result)
    if os.path.exists(hdock_out_path):
        print(f"Hdock output file already exists: {hdock_out_path}")
        hdock_score = read_hdock_scores(hdock_out_path)
        print(f"Hdock Score: {hdock_score}")
        return hdock_score

    with subprocess.Popen(
        # f"hdock {protein_path} {apta_path} -out {hdock_out_path}", # hdock aptamer_pdb/CCC.pdb protein_pdb_files/Seq_1.pdb -out hdock.out
        f'docker run --mount type=bind,source=C:\CodeProjects\github.com\JianJinglin\multiple-thread-hdock,target=/app win-hdock ./hdock/hdock {"/app/" + new_protein_path} {"/app/" + new_apta_path} -out {hdock_out_path}',
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True
    ) as process:
        while True:
            line = process.stdout.readline()
            if not line:
              break
            if 'seconds' in line:
              print(line, end='')
        process.wait()
        print(f"Finished hdock: {hdock_out_path}")
    hdock_score = read_hdock_scores(hdock_out_path)
    print(f"Hdock Score: {hdock_score}")
    return hdock_score


