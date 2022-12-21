# RNA-ERGY

This program contains various scripts for creating and training an objective function on a dataset of known (experimentally determined) 3D RNA structures. The objective function estimates the Gibbs free energy of a given structure, which can be used to evaluate and compare the stability of predicted RNA structures. In addition, this program includes scripts for generating scoring profiles, which plot the estimated Gibbs free energy as a function of interatomic distances. With these tools, you can optimize the 3D structure of RNA molecules and gain insights into their stability and function

## Installation 

```
git clone https://github.com/h-escoffier/RNA-ERGY
cd RNA-ERGY
pip install -r requirements.txt
```


## Usage

### training.py 

This script builds and trains an objective function using a dataset of 3D RNA structures that have been determined experimentally.

```
python3 trainig.py path_to_pdb_folder
```

### plot.py

This script generates scoring profiles, which plot the estimated Gibbs free energy as a function of interatomic distances. 

```
python3 plot.py 
or 
python3 plot.py path_to_energy_folder
```

### scoring.py 

This script uses the objective function to evaluate the accuracy of predicted structures from an RNA pdb file. 
```
python3 scoring.py path_to_pdb_file
or 
python3 scoring.py path_to_pdb_file path_to_energy_folder
```

