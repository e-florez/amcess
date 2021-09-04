from scipy.optimize import dual_annealing
from scipy.optimize import minimize, rosen, rosen_der
import numpy as np
nseed = np.random.seed(555)

def LennardJ(x):
	"""
	Lennard Jonnes en 1d
	e, profuncidad del pozo = 1
	sigma, distancia en la que el potencial es 0
	con los anteriores parámetros el mínimo está en -1
	"""
	#print("Lennard Jonnes 1D")
	E = 4.0*(x[0]**-12-x[0]**-6)
	boltzmann = np.exp(-E/273.15)
	return 4.0*(x[0]**-12 - x[0]**-6)

def LennardJv(x):
	"""
	Lennard Jonnes en 1d
	e, profuncidad del pozo = 1
	sigma, distancia en la que el potencial es 0
	con los anteriores parámetros el mínimo está en -1
	"""
	#print("Lennard Jonnes 1D")
	E = 4.0*(x**-12-x**-6)
	boltzmann = np.exp(-E/273.15)
	return 4.0*(x**-12 - x**-6)

def LJ3D(x):
	"""
	Lennard Jonnes para 2 átomos en 3D, donde uno de los 
	átomos está en el origen. En coordenadas Cartesianas
	"""
	print("Lennard Jones 3D",x)
	return 4.0*((x[0]**2+x[1]**2+x[2]**2)**-12 - (x[0]**2+x[1]**2+x[2]**2)**-6)

def Rf(x):
	"""
	Rastrigin function
	f(x) = An + Sum_i^n [x_i*x_i - Acos(2pix_i)]
	A = 10 y x e [-5.12,5.12]
	"""
	return np.sum(x*x - 10*np.cos(2*np.pi*x)) + 10*np.size(x)

def himmelblau(x):
	print("x  ",x) 
	return (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 -7)**2
#Raices
#f(3.0,2.0)=0.0,\quad
#f( − 2.805118 , 3.131312 )   = 0.0 
#f( − 3.779310 , − 3.283186 ) = 0.0 
#f( 3.584428 , − 1.848126 )   = 0.0

#bounds = [(-5,5),(-5,5)]
#bounds = [(-5,5),(-5,5),(-5,5)]

bounds = [(1,1.5)]
vx=[(1.0)]
#annealing = dual_annealing(func=LennardJ,bounds=bounds,local_search_options={'method': 'BFGS'}, maxiter=10000,initial_temp=100,restart_temp_ratio=0.6,seed=nseed,x0=vx,accept=-100) #,maxfun=100)
#annealing = dual_annealing(func=LennardJ,bounds=bounds,seed=nseed) #,maxfun=100)


#Funciona
#annealing = dual_annealing(func=LennardJ,bounds=bounds,local_search_options={'method': 'BFGS', 'jac' : 'cs', 'hess':'cs'})

#x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
#res = minimize(rosen, x0, args=0, method='dogleg', jac=rosen_der, hess='2-points',options={'gtol': 1e-6, 'disp': True})

#error couando uso jac, hes, cuando pongo más opciones en local_search_options

#maxfun : # veces que puedo llamar a func
#maxiter: # iteracciones por cada temperatura
#nfev   : # evaluaciones de la función objetivo
#njev   : # evaluaciones de la Jacobiana
#nhev   : # evaluaciones de la Hessiana

#local_search_options: Local minimización
# BFGS:   usa la primera derivada (local_search_options={'method': 'BFGS'})
# dogleg: usa gradiente y hessiana (local_search_options={'method': 'dogleg', 'jac':'BFGS', 'hess': 'HessianUpdateStrategy'})
# trust-ncg: Combina Newton metodología con Hessiana

lw = [-5.12] * 10 
up = [5.12] * 10
annealing = dual_annealing(LennardJ,list(zip(lw,up)), maxiter=1000, seed=nseed, local_search_options={'method': 'BFGS'})

for ix in annealing.x:
	print(LennardJv(ix),ix)
print("x: ",annealing.x[:],"f(x): ",annealing.fun,'Iteracciones ',annealing.nit)

#print(annealing)
#print(len(annealing.temperature))

#648     #Save initial temperature                                                        
#649     temperature_array = []
#650     temperature_array.append(initial_temp)
# En  vi /home/danian/.local/lib/python3.6/site-packages/scipy/optimize/_dual_annealing.py

