from src.base_molecule import Cluster
from src.search_configuration import SearchConfig


HF = {"atoms": [("H", 0, 0, 0), ("F", 0.917, 0, 0)]}
HF1 = {"atoms": [("H", 1, 8, 5), ("F", 0.918, 0.01, 10)]}

DHF = Cluster(HF, HF)

SC = SearchConfig(DHF)

SC.da()
