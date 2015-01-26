The aim of this project is to implement a **generator** which outputs **sequences** containing a **magic word**, and a **Gibbs sampler** which aims at finding the magic words in the previous sequences.

For this purpose, we generate N sequences si of length M, each containing a magic word of length w. Our Gibbs sampler should find R = {r1,...,rn} with ri the starting position of the magic word in the sequence i.

The **background** of a sequence is the union of every observation in a sequence which does not belong to the magic word.

A **magic word** is the union in a sequence of every observation having its index *i* between ri and ri + w - 1.

The experiment described in this report will focus on finding this motif described by a Dirichlet prior related to an alphabet. We have chosen to apply our Gibbs sampler to DNA sequences. Therefore, the alphabet K is K = {A, T, G, C}.
