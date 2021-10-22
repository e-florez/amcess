import src.m_GPyOpt
import src.math_cost_function

bounds = [(-6,6), (-6,6)]
opt = src.m_GPyOpt.solve_gaussian_processes(src.math_cost_function.Himmelblau, bounds, 5)
