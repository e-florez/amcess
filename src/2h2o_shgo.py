import sys
import numpy as np
sys.path.append('../')
from base_molecule import Cluster
from search_configuration import SearchConfig

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

################################################
# Configuración de la busqueda
###############################################
search_setting  = {
    "basis": 'sto-3g',
    "search_methodology": 2,
    "type of fragments": 0,
    "program_cost_function": 1,
    }

search_config = SearchConfig(h2o, search_setting)

search_config.run()
print(search_config._search)
