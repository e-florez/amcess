from amcess.base_molecule import Cluster, Molecule
from amcess.search_configuration import SearchConfig


HF = [("H", 0, 0, 0), ("F", 0.917, 0, 0)]
HF1 = {"atoms": [("H", 0, 0, 1), ("F", 0.918, 0.0, 1)]}
HF2 = {"atoms": [("H", 0, 0, 1), ("F", 0.918, 0.0, 1), ("He", 1, 2, 1)]}

DHF = Cluster(HF, HF1)

H2  = [("H", 0.0, 0.0, 0.0), ("H", 0.78, 0.0, 0.0)] 
H21 = [("H", 0.0, 0.0, 1.0), ("H", 0.78, 0.0, 1.0)]

#DHF = Molecule(H2)

SC = SearchConfig(DHF)

#SC.sphere_radius = 2.0
#SC.tolerance_contour_radius = 4.9

#SC.shgo()
SC.da()
#SC.bayesian()
