from amcess.base_molecule import Cluster
from amcess.search_configuration import SearchConfig


HF = [("H", 0, 0, 0), ("F", 0.917, 0, 0)]
HF1 = {"atoms": [("H", 1, 8, 5), ("F", 0.918, 0.01, 10)]}

DHF = Cluster(HF, HF)

SC = SearchConfig(DHF)

#SC.radius_contour = 2.0

#SC.shgo()
SC.da()
