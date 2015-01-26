import math
import operator
import functools
import itertools
import matplotlib.pyplot as plt
from numpy import random
from collections import Counter
from generator import generate_sample
from generator import generate_sequences

"""
This script generates a list of sequences, each containing a magic word (see generator.py)
and then estimates the most likely starting position of each of these motifs
"""

def estimate_position(alphabet, sequences, prev_positions, current_seq_idx, w, alpha_seq, alpha_motif):
    """
    Returns a vector of probabilities describing for each letter in the current sequence the probability
    for this position to be the starting position of the motif
    """
    p_positions = []
    N = len(sequences)
    M = len(sequences[0])
    alphabet_length = len(alphabet)

    for r_i in xrange(M-w):
        positions = list(prev_positions)
        positions[current_seq_idx] = r_i

        # Counter containing the nuber of times each character appears in the background
        backgrounds = [sequences[i][:positions[i]] + sequences[i][positions[i]+w:] for i in xrange(len(sequences))]
        c_background = Counter([item for sublist in backgrounds for item in sublist])

        # List of counters containing the number of times a character appears in the motifs for each column in a motif
        motifs = [sequences[i][positions[i]:positions[i]+w] for i in xrange(len(sequences))]
        c_motif = [Counter([motifs[i][j] for i in xrange(len(sequences))]) for j in xrange(w)]

        background_base = float(math.gamma(sum(alpha_seq))) / math.gamma(alphabet_length*(M-w) + sum(alpha_seq))
        motif_base = float(math.gamma(sum(alpha_motif))) / math.gamma(alphabet_length + sum(alpha_motif))

        # Compute probabilities
        prod = [float(math.gamma(c_background[alphabet[k]] + alpha_seq[k])) / math.gamma(alpha_seq[k]) for k in xrange(alphabet_length)]
        p_background = background_base * functools.reduce(operator.mul, prod, 1)

        for j in xrange(w):
            prod = [float(math.gamma(c_motif[j][alphabet[k]] + alpha_motif[k])) / math.gamma(alpha_motif[k]) for k in xrange(alphabet_length)]
            p_column = motif_base * functools.reduce(operator.mul, prod, 1)
            p_background *= p_column

        p_positions.append(p_background)

    return p_positions

def estimate_starting_positions(alphabet, sequences, alpha_seq, alpha_motif, w,
    iterations=200, minstep=100, steplength=10, plot=False):
    """
    Returns the estimated starting positions for the given sequences
    The results take into account only the estimations between minstep and iterations conditioned by steplength
    """
    positions = []
    N = len(sequences)
    M = len(sequences[0])

    # Initial positions randomly generated
    positions.append([random.randint(0, M-w+1) for i in xrange(N)])

    for i in xrange(iterations):
        tmp = []
        for j in xrange(N):
            p_pos = estimate_position(alphabet, sequences, positions[-1], j, w, alpha_seq, alpha_motif)

            # Normalization
            s = float(sum(p_pos))
            p_pos = [p / s for p in p_pos]

            # # Sampling
            pos = generate_sample(range(0, len(p_pos)), p_pos)

            tmp.append(pos)

        positions.append(tmp)

    # Plot convergence
    if plot:
        plt.figure(facecolor="white")
        plt.plot(positions)
        plt.title('Positions convergence for each sequence')
        plt.show()

    # We only select the iterations between minstep and iterations with a step of steplength
    positions = [positions[j] for j in xrange(minstep, iterations, steplength)]

    # For each sequence, the selected positions are the ones the most elected during the previous iterations
    return [Counter([positions[j][i] for j in xrange(len(positions))]).most_common(1)[0][0] for i in xrange(N)]


if __name__ == '__main__':
    alphabet = ['A', 'T', 'G', 'C']
    alpha_seq = [1,1,1,1]
    alpha_motif = [1,7,10,2]
    N = 5 # number of sequences
    M = 30 # sequence length
    w = 10 # motif length

    # Generate sequences and positions
    sequences, positions = generate_sequences(alphabet, alpha_seq, alpha_motif, N, M, w)

    # Estimates the starting positions
    iterations = 200
    minstep = 100
    steplength = 10
    estimated_positions = estimate_starting_positions(alphabet, sequences, alpha_seq, alpha_motif, w,
        iterations=iterations, minstep=minstep, steplength=steplength, plot=True)

    print estimated_positions

    print 'Results for N={}, M={}, w={}, nb_iter={}, min_step={}, stepsize={}, alpha_motif={}'.format(
        N, M, w, iterations, minstep, steplength, alpha_motif)
    for (pos, estimated_pos) in zip(positions, estimated_positions):
        print '{} -> {} ({})'.format(pos, estimated_pos, 'OK' if pos == estimated_pos else 'FALSE')
