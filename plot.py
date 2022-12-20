import numpy as np
import matplotlib.pyplot as plt
import os
import argparse


def interaction_profile_plot(file):
    """
    Create a plot representing an interaction profile, pairwise scores as a function of distance (Å), taking a file
    as input
    :param file: File containing the pairwise score. Each line of a file represents an interval of 1Å. (e.g. line 1
    represents the interval [0, 1] Å.)
    """
    x_values = [i for i in range(20)]
    y_values = []
    with open('energy/' + file, 'r') as pb_file:
        content = [x.strip('\n') for x in pb_file.readlines()]
    for elm in content:
        if elm == 'None':
            y_values.append(np.nan)
        else:
            y_values.append(float(elm))
    fig, ax = plt.subplots()
    ax.plot(x_values, y_values, linewidth=2.0)
    plt.title('Interaction profile')
    plt.legend([file])
    ax.set_xlabel('Distance (Å)')
    ax.set_ylabel('Pairwise score')
    ax.set(xlim=(0, 20))
    plt.savefig('plot/' + str(file) + '.png')


def main(folder='Energy'):
    """
    Create graphs of the interaction profiles for each base pair and save them in a 'plot' folder.
    :param folder: Folder containing 10 file for each base pair containing for each the pairwise scores.
    """
    try:  # Create the 'plot' folder if it does not exist
        os.mkdir('plot')
    except OSError:
        pass
    for bp in os.listdir(folder):
        interaction_profile_plot(bp)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('energy_folder', nargs='?', help="Give the path of the folder containing the energy file for the 10 base pairs. By default the name of this file is 'Energy'.")
    args = parser.parse_args()
    if args.energy_folder:
        main(args.energy_folder)
    else:
        main()
