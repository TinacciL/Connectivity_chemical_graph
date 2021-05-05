# Connectivity in chemical graph
This is a Python 3.8 script that permits to encode molecules as molecular graph and compare the connectivity of two molecules.

## Installation

This program used a python3 interface, to run this code you must install on your machine this list of packages:

* ```networkx```
* ```ase```

## Script Structure

0 - Import the required packages

```python
import networkx as nx
from ase import io, neighborlist, atoms
```

1 - Function that encode the molecule from the Euclidean space (.xyz) to the corrisponding graph structure, and you will call it in the main of the code.

```python
def FromXYZtoGraph(input_file):
    atoms = ['H','He','Li','C','N','O','F','Na','Si','P','S','Cl']
    atomic_numb = [1,2,3,6,7,8,9,11,14,15,16,17]
    mol = io.read(input_file)
    #compute neighbor of the atoms in xyz format
    cutOff = neighborlist.natural_cutoffs(mol)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
    neighborList.update(mol)
    #compure adjacency matrix and atoms list
    adj_matrix = neighborList.get_connectivity_matrix(sparse=False)
    Natom_list = mol.get_atomic_numbers()
    atoms_list = []
    for i,item in enumerate(Natom_list):
        for k in range(len(atomic_numb)):
            if item == atomic_numb[k]:
                atoms_list.append(atoms[k]) 
    #convert in networkx-molecules graph
    G=nx.from_numpy_matrix(adj_matrix)
    for i,item in enumerate(atoms_list):
        tmp_attr = {'atom': item}
        G.nodes[i].update(tmp_attr.copy())
    return(G)
```
2 - Code main: computing the isomorphism between two molecular graph object.

   2.1 - Encode the molecoles from .xyz to graph molecule objects, modify the path (example "pathToMolecules/mol_0.xyz") where your .xyz is locate.
```python
mol_0 = FromXYZtoGraph(pathToMolecules/mol_0.xyz)
mol_1 = FromXYZtoGraph(pathToMolecules/mol_1.xyz)
```
    2.2 - Run the function to control the isomorphism between the two molecules
```python
prop = 'atom'
nm = nx.algorithms.isomorphism.categorical_node_match(prop,prop)

if nx.is_isomorphic(mol0,mol_1,node_match=nm):
    print('Are isomorphic!')
else:
    print('Are NOT isomorphic!')
```

## Associated publication
aaa

## Acknowledgments
This project has received funding within the European Union’s Horizon 2020 research and innovation programme from the Marie Sklodowska-Curie for the project ”Astro-Chemical Origins” (ACO), grant agreement No 811312.
