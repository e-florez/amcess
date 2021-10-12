import sys

sys.path.append("../")

from src import base_molecule, search_configuration


#################################################
# Configuración del sistema de trabajo
#################################################
water1 = {
    "atoms": [
        ("H", 1.4891244, 0.0, 1.0332262),
        ("O", 0.0, 0.0, -0.1302053),
        ("H", -1.4891244, 0.0, 1.0332262),
    ]
}

h2o = base_molecule.Cluster(water1, water1)
# Configuración de busqueda
################################################
# Configuración de la busqueda
###############################################


search_config = search_configuration.SearchConfig(h2o)

search_config.da()
print(search_config._search)
