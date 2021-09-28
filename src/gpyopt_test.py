import GPyOpt

def himmelblau(x0) -> float:
    """ Himmelblau's function is a multi-modal function, used to test the 
    performance of optimization algorithms, 
    
            f(x, y) = (x**2 + y - 11)**2 + (x + y**2 - 7)**2

    Parameters
    ----------
    x : float
        [description]
    y : float
        [description]

    Returns
    -------
    float
        [description]
    """
    x = x0[0,0]
    y = x0[0,1]
    return ((x**2 + y - 11)**2 + (x + y**2 -7)**2)

initer = 20
minval = -6
maxval = 6
bounds = [{'name': 'x', 'type': 'continuous', 'domain': (minval, maxval)},
          {'name': 'y', 'type': 'continuous', 'domain': (minval, maxval)}]
xi = 0.001

objective = GPyOpt.core.task.SingleObjective(himmelblau)
space = GPyOpt.Design_space(space = bounds)
model = GPyOpt.models.GPModel(optimize_restarts=5,verbose=False)
acquisition_optimizer = GPyOpt.optimization.AcquisitionOptimizer(space)
initial_design = GPyOpt.experiment_design.initial_design('latin', space, initer)
# acquisition = jitter_integrated_EI(model, space, optimizer=aquisition_optimizer, par_a=1, par_b=10, num_samples=200)
acquisition=GPyOpt.acquisitions.AcquisitionEI(model, space, acquisition_optimizer,jitter=xi)
evaluator = GPyOpt.core.evaluators.Sequential(acquisition)
opt_himmelblau = GPyOpt.methods.ModularBayesianOptimization(model, space, objective, acquisition, evaluator, initial_design)

maxiter = 50
# opt_himmelblau = GPyOpt.methods.BayesianOptimization(himmelblau, domain=bounds, initial_design_numdata=initer, initial_design_type='latin')
# opt_himmelblau.run_optimization(max_iter = maxiter, save_models_parameters=True, report_file='report.log', evaluations_file='evals.log', models_file='model.log')
opt_himmelblau.run_optimization(max_iter = maxiter, save_models_parameters=True, evaluations_file='evals.log', models_file='model.log')
opt_himmelblau.plot_acquisition()
opt_himmelblau.plot_convergence()
print(opt_himmelblau.x_opt, opt_himmelblau.fx_opt)