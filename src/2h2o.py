import sys
import numpy as np
sys.path.append('../')
from base_molecule import Cluster
#from m_dual_annealing import solve_dual_annealing
from search_configuration import SearchConfig

#Sistema de estudio
#hydrogen2 = {  # bond distance = 74 pm
#    "coordinates": [
#        ("H", 1, 1, 0),
#        ("H", 1.74, 1, 0),
#    ],
#    "charge": 0,
#    "multiplicity": 1,
#}
#################################################
# Configuración del sistema de trabajo
#################################################
water1 = {
    "coordinates": [
        ("H", 1.4891244, 0., 1.0332262),
        ("O", 0., 0., -0.1302053),
        ("H", -1.4891244, 0., 1.0332262),
    ]
}
water2 = {
    "coordinates": [
        ("H", 2.4891244, 1., 2.0332262),
        ("O", 1., 1., 0.8697947),
        ("H", -0.4891244, 1., 2.0332262),
    ],
}

h2o = Cluster(water1, water2)

#x0 = np.asarray(h2o.center_of_mass)
#print(h2o.center_of_mass)

x0 = [0.00000000e+00,  0.00000000e+00, -9.46852623e-06, 1., 1., 0.99999053]

#Configuración de busqueda
a = 1.0
################################################
# Configuración de la busqueda
###############################################
search_setting  = {
    "bounds": [(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a)],
    "basis": 'sto-3g',
    "search_methodology": 1,
    "type of fragments": 0,
    "program_cost_function": 1
    }

search_config = SearchConfig(h2o, search_setting)

#Empieza la busqueda
#parameters = {"NT":3,
#            "no_local_search":True}

#search_config.run(no_local_search=False, args=('H', 'O', 'H', 'H', 'O', 'H', main_axes, water_dimer),
search_config.run(no_local_search=False,
                  x0 = x0)
exit()
exit()
#search_config.run(NT=3, no_local_search=True)
print(search_config._search)



print("search_config \n",search_config)


print(f"\t{h2.total_atoms}")
print(f" rot = 0 deg")
for i in range(h2.total_atoms):
    print("\t".join(h2.write_molecule[i]))


for i in range(8):
    ang = 45 * (i + 1)
    h2.coordinates = h2.rotate(x=0, y=0, z=45)
    h2.coordinates = h2.translate(x=0.4, y=0, z=0)

    print(f"\t{h2.total_atoms}")
    print(f" rot = {ang} deg")
    for i in range(h2.total_atoms):
        print("\t".join(h2.write_molecule[i]))
