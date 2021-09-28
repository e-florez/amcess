import sys
import numpy as np
sys.path.append('../')
from base_molecule import Cluster
#from m_dual_annealing import solve_dual_annealing
from search_configuration import SearchConfig

#Sistema de estudio
hydrogen2 = {  # bond distance = 74 pm
    "coordinates": [
        ("H", 1, 1, 0),
        ("O", 0,0,0),
        ("H", 2.74, 1, 1),
    ],
    "charge": 0,
    "multiplicity": 1,
}

h2 = Cluster(hydrogen2)

x0 = [1, 1, 0, 0, 0, 0, 2.74, 1, 1]
print(x0)
#Configuración de busqueda
a = 2.5
################################################
# Configuración de la busqueda
###############################################
search_setting  = {
    "bounds": [(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a)],
    "basis": 'sto-3g', #pyscf
    "search_methodology": 0,
    "type of fragments": 0,
    "program_cost_function": 0, #0: pyscf, 1: ps4, 2: horton, 3: gaussian, ...
    }
search_config = SearchConfig(h2, search_setting)
print(search_config._search_name)

#Empieza la busqueda

search_config.run(no_local_search=False,
                  x0 = x0)
print(search_config._search)
exit()



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
