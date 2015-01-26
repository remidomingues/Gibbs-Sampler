from math import ceil, trunc
import matplotlib.pyplot as plt
from itertools import cycle
from generator import generate_sequences
from matplotlib.font_manager import FontProperties
from gibbs_sampler import estimate_starting_positions

"""
This scripts plot the convergence of the positions for each sequence (Markov chain generated)
and the accuracy statistics for various parameters
"""

def plot_statistics(statistics, x_label='X', y_label='Y', fct_label='', fig_param_label=''):
    """
    statistics: [(fct_param, x, fig_param, y)]
    """
    fig_params = sorted(set([fig_param for (fct_param, x, fig_param, y) in statistics]))
    fct_params = sorted(set([fct_param for (fct_param, x, fig_param, y) in statistics]))
    x_params = sorted(set([x for (fct_param, x, fig_param, y) in statistics]))
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')

    # We regroup statistics by fig_param value
    fig_params_stats = [[(fct_param, x, y) for (fct_param, x, fig, y) in statistics if fig == fig_param] for fig_param in fig_params]

    # One figure per fig_param
    for i in xrange(len(fig_params)):
        # We regroup statistics for a given fig_param by fct_params value
        fct_params_stats = [[(x, y) for (fct, x, y) in fig_params_stats[i] if fct == fct_param] for fct_param in fct_params]

        functions = []
        plt.figure(facecolor="white")
        ax = plt.subplot(111)
        ax.axis([x_params[0], x_params[-1], 0, 1])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.8])
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title('{} ({}={})'.format(y_label, fig_param_label, fig_params[i]))

        # One function per fct_params
        for j, color in zip(xrange(len(fct_params)), colors):
            # X is x_params
            _x = [x for (x, y) in fct_params_stats[j]]
            # Y is classification rate
            rates = [y for (x, y) in fct_params_stats[j]]
            fct, = ax.plot(_x, rates, color, label='{}={}'.format(fct_label, fct_params[j]), lw=2)
            ax.plot(_x, rates)
            functions.append(fct)

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
                  fancybox=True, shadow=True, ncol=5, handles=functions)
        plt.show()


if __name__ == '__main__':
    statistics = []
    alphabet = ['A', 'T', 'G', 'C']
    alpha_seq = [1,1,1,1]
    alpha_motif = [1,7,20,2]
    N = 5 # number of sequences
    M = 30 # sequence length
    w = 10 # motif length

    # Run benchmarks for results estimation
    # Generate sequences and positions
    sequences, positions = generate_sequences(alphabet, alpha_seq, alpha_motif, N, M, w)

    # # Run benchmarks
    # for iterations in [10,30]:
    #     print 'iterations=' + str(iterations)
    #     for minstep in range(0, iterations/2, iterations / 5):
    #         for stepsize in range(1,10,2):
    #             # Estimate the starting positions
    #             estimated_positions = estimate_starting_positions(alphabeit, sequences, alpha_seq, alpha_motif, w,
    #                 iterations=iterations, minstep=minstep, steplength=stepsize)

    #             # Compute accuracy
    #             accuracy = sum([1 if pos == estimated_pos else 0
    #                 for pos, estimated_pos in zip(positions, estimated_positions)]) / float(len(positions))
    #             print accuracy
    #             statistics.append((stepsize, minstep, iterations, accuracy))

    # # Plot statistics
    # plot_statistics(statistics, x_label='min_step', y_label='Accuracy',
    #     fct_label='step_size', fig_param_label='Number of iterations')


    # Run benchmarks for problem complexity
    iterations = 100
    minstep = 50
    stepsize = 5
    nbruns = 5

    for N in range(1, 10, 1):
        print 'N=' + str(N)
        for M in [20]:
            for w in range(2, M, int(M*0.2)):
                accuracy = 0
                for i in xrange(nbruns):
                    # Generate sequences and positions
                    sequences, positions = generate_sequences(alphabet, alpha_seq, alpha_motif, N, M, w)

                    # Estimate the starting positions
                    estimated_positions = estimate_starting_positions(alphabet, sequences, alpha_seq, alpha_motif, w,
                        iterations=iterations, minstep=minstep, steplength=stepsize)

                    # Compute accuracy
                    accuracy += sum([1 if pos == estimated_pos else 0
                        for pos, estimated_pos in zip(positions, estimated_positions)]) / float(len(positions))

                accuracy /= float(nbruns)
                print accuracy
                statistics.append((w, N, M, accuracy))

    # Plot statistics
    plot_statistics(statistics, x_label='Number of sequences', y_label='Accuracy',
        fct_label='w', fig_param_label='Sequence length M')
