from amcess.base_molecule import Cluster, Molecule
from amcess.search_engine import SearchConfig


HF = [("H", 0, 0, 0), ("F", 0.917, 0, 0)]
HF1 = {"atoms": [("H", 0, 0, 1), ("F", 0.918, 0.0, 1)]}
w = [("O", 0, 0, 0),("H", 0.58708, 0.75754, 0),("H", -0.58708, 0.75754, 0)]


#DHF = Cluster(HF, HF1)
w4 = Cluster(w,w,w,w)

w4 = w4.initialize_cluster()

#w4.sphere_center = (0,0,0)
#w4.sphere_radius = 2.0

SC = SearchConfig(w4)
SC._methodology = "MP2"
#SC._basis_set = "ccpvdz"
# print(SC.search_type)
#SC.search_type = "Bayesian"
#SC.search_type = "SHGO"
#SC.search_type = "dual_annealing"
# SC.search_type = "Edison"
# SC.sphere_radius = 2.0
# SC.tolerance_contour_radius = 4.9
#print("centro ",SC._sphere_center)
#print("radio  ",SC._sphere_radius)
#exit()
#SC.run(gp_params={"initer": 3, "maxiter": 3})  # Bayesian
#SC.run(initer=10, maxiter=10)  # Bayesian
#SC.run(sampling_method="sobol", n=2)
#SC.run(maxfun=10,maxiter=10)  # !da
SC.run(T0=1000, dT=0.2, nT=20, maxCycle=100)  #
#SC.run()
