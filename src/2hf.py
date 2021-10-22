import sys

sys.path.append('../data')

from base_molecule import Cluster, Molecule
from search_configuration import SearchConfig


HF={"atoms": [("H", 1, 3, 4), ("F", 0.917, 0, 0)]}
HF1={"atoms": [("H", 1, 8, 5), ("F", 0.918, 0.01, 10)]}

DHF = Cluster(HF1, HF1)

SC = SearchConfig(DHF)

