from Bio import Phylo
trees = Phylo.read('final_super_alignment.phy_phyml_tree.txt', 'newick')
print(trees)
Phylo.draw(trees)
