from training import Cleaner
from shared import bp_attribution
import sys


def main(pdb_file):
    """
    Calculation of the Gibbs free energy of an RNA input (in .pdb format) using the objective function trained
    in the script 'training.py'
    :param pdb_file:
    :return:
    """
    c3_list = Cleaner(pdb_file).run()
    gibbs_free_energy = bp_attribution(c3_list, 'scoring', gibbs_free_energy=0)
    print('Gibbs Free Energy: ', gibbs_free_energy)


if __name__ == '__main__':
    try:
        main(sys.argv[1])
    except IndexError:
        print('Give the path of a PdB file.')
