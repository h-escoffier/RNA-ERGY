from training import Cleaner
from shared import bp_attribution
import argparse


def main(pdb_file, energy='Energy'):
    """
    Calculation of the Gibbs free energy of an RNA input (in .pdb format) using the objective function trained
    in the script 'training.py'
    :param energy: Path to a folder 'Energy' containing 10 files (one per base pair) of 20 lines each.
    Each line contains the pseudo energy of an interval of 1 between 0 and 20.
    :param pdb_file:
    :return:
    """
    c3_list = Cleaner(pdb_file).run()
    gibbs_free_energy = bp_attribution(c3_list, 'scoring', gibbs_free_energy=0, energy=energy)
    print('Gibbs Free Energy: ', gibbs_free_energy)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_file', help="Give the path of a PdB file")
    parser.add_argument('energy', nargs='?', help="Give the path to the folder Energy")
    args = parser.parse_args()
    if args.energy:
        main(args.pdb_file, args.energy)
    else:
        main(args.pdb_file)

