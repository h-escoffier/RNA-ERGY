from shared import Cleaner, bp_attribution
import numpy as np
import os
import sys


def create_freq_list():
    """
    Create the list in which the frequencies for each bp will be stored.
    The list is of the form [[[AB], [0,0...0]], [[CD], [0,0...0]] ... [[YZ], [0,0...0]]].
    The last element of the list 'XX' contains the reference frequency.
    :return: The list for each frequency for each bp (and the reference) set to 0
    """
    frequency = []
    for bp in ['AA', 'AU', 'AC', 'AG', 'UU', 'CU', 'GU', 'CC', 'CG', 'GG', 'XX']:
        frequency.append([bp, 20 * [0]])
    return frequency


def obs_ref_frequency(frequency):
    """
    Convert the number of elements in each interval to the observed frequency for each interval.
    The last element in the list 'XX' is the reference frequency
    """
    for pb in frequency:
        if sum(pb[1]) != 0:
            pb[1] = [nr / sum(pb[1]) for nr in pb[1]]
    return frequency


def pseudo_energy(score, frequency):
    """
    Calculates the score (pseudo-energy) for each interval of each base pair. 
    """
    for pb in frequency:
        pe = []
        index = 0
        for fobs in pb[1]:
            if fobs != 0 and frequency[10][1][index] != 0:
                pseudo_score = -np.log(fobs / frequency[10][1][index])
                if pseudo_score > 10:
                    pe.append(10)
                else:
                    pe.append(pseudo_score)
            else:
                pe.append(10)
            index += 1
        score.append([pb[0], pe])
    score.pop()
    return score


def file_generator(score):
    """
    Create 10 files (one per base pair) of 20 lines each. Each line contains the pseudo energy of an interval of 1 between 0 and 20.
    :param score: A list of 10 lists of the form : [[[AA], [0,0...0]], [[AC], [0,0...0]] ... [[UU], [0,0...0]]]
    An element as the form: [[AA][x,y, ... , z]] with x,y,z respectively the pseudo-energy of the intervals [0,1], [1,2] and [19,20].
    """
    try:  # Create the 'energy' folder if it does not exist
        os.mkdir('Energy')
    except OSError:
        pass
    for pb in score:
        with open('energy/' + pb[0], 'w+') as pb_file:
            content = [str(x) for x in pb[1]]
            pb_file.write("\n".join(content))


def main(pdb_folder):
    """
    Calculate the inter-atomic distances of each bp from a given set of '.pdb' files in order to classify them in
    intervals from 0 to 20, then convert these distances into pseudo-energy and finally generate 10 files (1 for each
    bp) of 20 lines (each line corresponding to an interval)
    :param pdb_folder: Folder containing one or more '.pdb' files
    """
    frequency = create_freq_list()
    score = []
    for pdb_file in os.listdir(pdb_folder):  # Browse all the '.pdb' files in the folder
        c3_list = Cleaner(pdb_folder + '/' + pdb_file).run()
        frequency = bp_attribution(c3_list, 'training', frequency=frequency)
    frequency = obs_ref_frequency(frequency)
    score = pseudo_energy(score, frequency)
    file_generator(score)


if __name__ == '__main__':
    try:
        main(sys.argv[1])
    except IndexError:
        print('Give the path of the folder containing the PdB Files as an argument.')
