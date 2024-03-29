#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the
#   AMCESS Project (https://gitlab.com/ADanianZE/amcess).
# Copyright (c) 2021, Edison Florez
# License: GPLv3
#   Full Text: https://gitlab.com/ADanianZE/amcess/-/blob/main/LICENSE

# =============================================================================
# DOCS
# =============================================================================

"""
AMCESS
Atomic and Molecular Cluster Energy Surface Sampler.
"""

# =============================================================================
# META
# =============================================================================

__name__ = "amcess"
__version__ = "0.1.2a20"


# =============================================================================
# IMPORTS
# =============================================================================

from .base_molecule import Molecule, Cluster  # noqa
from .search_engine import SearchConfig as engine  # noqa
