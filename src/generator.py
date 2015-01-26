from numpy import random
from collections import Counter
"""
This script generates N sequences of length M, each containing a magic word of length w
The letters in the sequences and magic words are obtained from a categorical distribution
from Dir(alpha), with alpha a vector of reals of the size of the alphabet
"""

def generate_sample(alphabet, categ):
    """
    Returns a sample from the given categorical distribution
    """
    r = random.random()
    idx = 0
    s = categ[idx]

    while s < r:
        idx += 1
        if idx == len(categ):
            return alphabet[len(alphabet)-1]
        s += categ[idx]

    return alphabet[idx]

def generate_sequences(alphabet=['A', 'T', 'G', 'C'], alpha_seq=[1,1,1,1], alpha_motif=[1,7,12,2], N=5, M=10, w=3):
    """
    Returns a list of sequences and the list of magic words starting position
    """
    # Sequence generation
    categ_seq = random.dirichlet(alpha_seq)
    sequences = [[generate_sample(alphabet, categ_seq) for j in xrange(M)] for i in xrange(N)]

    # Magic words generation
    categs_motif = random.dirichlet(alpha_motif, w)
    positions = [random.randint(0, M-w+1) for i in xrange(N)]
    sequences = [[generate_sample(alphabet, categs_motif[idx-pos]) if idx >= pos and idx < pos+w else seq[idx] for idx in xrange(len(seq))]
        for seq, pos in zip(sequences, positions)]

    return sequences, positions


if __name__ == '__main__':
    # generate_sample test
    categ = [0.1, 0.6, 0.1, 0.2]
    values = ['A', 'B', 'C', 'D']
    c = Counter([generate_sample(values, categ) for i in xrange(1000)])
    print 'test sampling: ' + str(c)

    # generate_sequences test
    sequences, positions = generate_sequences()
    print 'sequences: ' + str(sequences)
    print 'positions: ' + str(positions)
