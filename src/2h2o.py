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
    "atoms": [
        ("H", 1.4891244, 0., 1.0332262),
        ("O", 0., 0., -0.1302053),
        ("H", -1.4891244, 0., 1.0332262),
    ]
}
water2 = {
    "atoms": [
        ("H", 2.4891244, 1., 2.0332262),
        ("O", 1., 1., 0.8697947),
        ("H", -0.4891244, 1., 2.0332262),
    ],
}

h2o = Cluster(water1, water2)

#Configuración de busqueda
a = 0.5
################################################
# Configuración de la busqueda
###############################################
search_setting  = {
    "bounds": [(-a, a),(-a, a),(-a, a),(-a, a),(-a, a),(-a, a)],
    "basis": 'sto-3g',
    "search_methodology": 1,
    "type of fragments": 0,
    "program_cost_function": 1,
    }

search_config = SearchConfig(h2o, search_setting)

search_config.run(x0 = [0,0,0,0,0,0])
print(search_config._search)
