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
water_dimer = {
    "coordinates": [
        ("H", 1.4891244, 0., 1.0332262),
        ("O", 0., 0., -0.1302053),
        ("H", -1.4891244, 0., 1.0332262),
        ("H", 2.4891244, 1., 2.0332262),
        ("O", 1., 1., -2.1302053),
        ("H", -2.4891244, 1., 2.0332262),
    ],
    "charge": 0,
    "multiplicity": 1,
}

h2o = Cluster(water_dimer)

x0 = np.asarray(h2o.coordinates)

centermass = h2o.center_of_mass + h2o.center_of_mass
main_axes = h2o.principal_axis + h2o.principal_axis
#Configuración de busqueda
a = 1.0
################################################
# Configuración de la busqueda
###############################################
search_setting  = {
    "bounds": [(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),
    (-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),
    (-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a)],
    "basis": 'sto-3g',
    "search_methodology": 0,
    "type of fragments": 1,
    "program_cost_function": 0
    }

search_config = SearchConfig(h2o, search_setting)
exit()
#print(search_config._bounds)
#print(search_config._basis)
search_config.search_name()
print(search_config._search_name)

#Empieza la busqueda
#parameters = {"NT":3,
#            "no_local_search":True}
x0 = np.zeros((6),dtype=float)
#x0 = np.asarray(centermass) #pys tiene conflictos exit cuando moleculas superpuesta
#x0 = [0,0,0,1,1,1]
x0 = [1.4891244, 0., 1.0332262, 0., 0., -0.1302053, -1.4891244, 0., 1.0332262,
    2.4891244, 1., 2.0332262, 1., 1., 1.-0.1302053, 1.-1.4891244, 1., 2.0332262,]

#search_config.run(no_local_search=False, args=('H', 'O', 'H', 'H', 'O', 'H', main_axes, water_dimer),
search_config.run(no_local_search=False,
                  x0 = x0)
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
