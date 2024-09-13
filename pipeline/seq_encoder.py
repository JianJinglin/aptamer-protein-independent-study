pair_to_char = {
    'AA': 'Q', 'AU': 'W', 'AC': 'E', 'AG': 'R',
    'UA': 'T', 'UU': 'Y', 'UC': 'U', 'UG': 'I',
    'CA': 'O', 'CU': 'P', 'CC': 'A', 'CG': 'S',
    'GA': 'D', 'GU': 'F', 'GC': 'G', 'GG': 'H',
    'AN': 'J', 'UN': 'K', 'CN': 'L', 'GN': 'M',
    'NA': 'Z', 'NU': 'V', 'NC': 'B', 'NG': 'N',
}

char_to_pair = {v: k for k, v in pair_to_char.items()}

def encode_rna(seq):
    if len(seq) % 2 != 0:
        seq += 'N'
    encoded = ''.join(pair_to_char[seq[i:i+2]] for i in range(0, len(seq), 2))
    return encoded

def decode_rna(encoded):
    decoded = ''.join(char_to_pair[char] for char in encoded)
    if decoded.endswith('N'):
        decoded = decoded[:-1]
    return decoded