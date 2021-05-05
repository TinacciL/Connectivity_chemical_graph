import networkx as nx
from ase import io, neighborlist, atoms

#function that encode the molecules from xyz file to graph object
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

#encode the molecules from xyz file to graph object
mol_0 = FromXYZtoGraph(pathToMolecules/mol_0.xyz)
mol_1 = FromXYZtoGraph(pathToMolecules/mol_1.xyz)

#test the isomorphism
prop = 'atom'
nm = nx.algorithms.isomorphism.categorical_node_match(prop,prop)

if nx.is_isomorphic(mol0,mol_1,node_match=nm):
    print('Are isomorphic!')
else:
    print('Are NOT isomorphic!')
