Dual Annealing
--------------

In addition, the Dual Annealing method, implemented in the `Scipy Library <https://docs.scipy.org/doc/scipy/index.html#>`_, has also been included into the AMCESS package. This method follows the Generalized Simulated Annealing, i.e., a generalization of the SA :cite:`tsallis1988possible,tsallis1996generalized,xiang1997generalized,xiang2000efficiency,xiang2013generalized,mullen2014continuous`. This technique uses the Tsallis-Stariolo form of the Cauchy-Lorentz visiting distribution, which avoids the possibility of getting trapped in local minima so the global minimum can be found. Moreover, a generalized Metropolis algorithm is used for the acceptance probability, which also provides the possibility of finding multiple local minima, if there are any.
