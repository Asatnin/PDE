import numpy as np
from scipy.optimize import curve_fit


def exponential_fit(x, a, b, c):
    return a*np.exp(-b*x) + c


x = np.array([293.0, 1278.0, 1528.0, 1677.0])
y = np.array([2.0 * 0.001, 5.0 * 0.001, 7.8 * 0.001, 0.01])
fitting_parameters, covariance = curve_fit(exponential_fit, x, y)
a, b, c = fitting_parameters


def next_k_p(t):
    return exponential_fit(t, a, b, c)