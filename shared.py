import math
import numpy as np


class Cleaner:
    # This class takes as input a '.pdb' file and outputs a list containing only the 'c3' atoms in that file.
    def __init__(self, path_pdb_file):
        self.path_pdb_file = path_pdb_file
        self.atom_list = []
        self.c3_list = []

    # 'Return' a List containing all atoms
    def atom_extractor(self):
        pdb_file = open(self.path_pdb_file, 'r')
        for line in pdb_file:
            if line[0:4] == 'ATOM':
                cols = line.strip('\n')
                # Handles the particular structure of the '.pdb' file
                atom = [cols[0:6], cols[6:11], cols[12:16], cols[16:17], cols[17:20], cols[21:22], cols[22:26],
                        cols[26:27], cols[30:38], cols[38:46], cols[46:54], cols[54:60], cols[60:66], cols[76:78],
                        cols[78:80]]
                elm = [i.replace(' ', '') for i in atom]
                self.atom_list.append(elm)
        pdb_file.close()

    # Select only c3' atoms of the atom_list
    def c3_extractor(self):
        for elm in self.atom_list:
            if elm[2] == "C3'":
                self.c3_list.append(elm)

    def run(self):
        self.atom_extractor()
        self.c3_extractor()
        return self.c3_list


def distance_calculator(res1, res2):
    """
    Calculates the distance between two residuals
    :param res1: Residual 1
    :param res2: Residual 2
    :return: A float representing the distance between the two residuals taken as input
    """
    res1 = list(np.float_(res1))
    res2 = list(np.float_(res2))
    distance = ((res1[0] - res2[0]) ** 2 + (res1[1] - res2[1]) ** 2 + (res1[2] - res2[2]) ** 2) ** (1 / 2)
    return distance


def bp_attribution(c3_list, objective, frequency=None, gibbs_free_energy=None):
    """
    From a list of bases one creates pairs according to pre-established conditions and either adds them according to
    the distance and the base pair to a list or calculates the Gibbs free energy
    :param c3_list: List where each element is a line of a C3' atom extracted from a pdb file and preprocessed by the Cleaner class.
    :param objective: 'training' or 'scoring' depending on whether the function is used to train the objective function (training)
    or to estimate the Gibbs free energy of an RNA conformation (scoring).
    :param frequency: Only for 'training'
    The 'frequency' list is of the form [[[AA], [0,0...0]], [[AC], [0,0...0]] ... [[XX], [0,0...0]]] with the last element of the list 'XX' contains the reference frequency.
    In this part the corresponding distance interval is incremented by 1 for each base pair.
    :param gibbs_free_energy: Only for 'scoring'
    :return:
    """
    for res1 in c3_list:
        compare = False
        for res2 in c3_list:
            if compare:  # Allows a base pair to be processed only once
                if res1[5] == res2[5] and (int(res1[6]) + 4) <= int(res2[6]):  # Conditions: Only 'intra-chain' distances are considered & Only consider residues separated by at least 3 positions on the sequence
                    bp = ''.join(sorted((res1[4] + res2[4])))
                    dist = distance_calculator(res1[8:11], res2[8:11])
                    if objective == 'training':
                        if 0 <= dist <= 20:
                            for bp_tot in frequency:
                                if bp == bp_tot[0]:
                                    int_dist = math.floor(dist)
                                    bp_tot[1][int_dist] += 1
                                    frequency[10][1][int_dist] += 1
                    else:  # type == 'scoring'
                        if 0.5 <= dist <= 19.5:
                            gibbs_free_energy += linear_interpolation(bp, dist)
            if res1 == res2:
                compare = True
    if objective == 'training':
        return frequency
    else:
        return gibbs_free_energy


def linear_interpolation(bp, dist):
    """
    A score value is calculated using linear interpolation.
    Each line is taken as the average of the interval. (e.g. line 1 corresponding to the interval [0, 1] is
    considered to be the value 0.5 in the interpolation.)
    Thus, only values between 0.5 and 19.5 are taken into account.
    Note this program is only used in the 'scoring.py' script but is in this file to avoid an error due to a circular import.
    :param bp: A base pair
    :param dist: Distance between the two bases
    :return: Energy associated with base pair and distance
    """
    with open('Energy/' + bp, 'r') as energy:
        content = energy.readlines()
    energy_before = float(content[math.floor(dist) - 1])
    energy_after = float(content[math.floor(dist)])
    energy = energy_before + (dist - (math.floor(dist) - 0.5)) * (energy_after - energy_before)  # Linear Interpolation Formula
    return energy
