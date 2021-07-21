import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

font_params = {
    'font.size': 10,
    'axes.labelsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
}

plt.rcParams.update(**font_params)


def PolyCoefficients(x, coeffs):
    """ Returns a polynomial for ``x`` values for the ``coeffs`` provided.

    The coefficients must be in ascending order (``x**0`` to ``x**o``).
    """
    o = len(coeffs)
    print(f'# This is a polynomial of order {ord}.')
    y = 0
    for i in range(o):
        y += coeffs[i] * x**i
    return y


def monoExp(x, m, t, b):
    return m * np.exp(-t * x) + b


def powerlaw(x, m, t):
    return m * np.power(x, t)


def linlaw(x, m, t, ):
    return t * x + m
# def powerlaw(x, m, t):
#     return m * x**(-t)


def fit_power(x, y, p0, fit_func):
    # p0 = (0, 0)  # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(fit_func, x, y, p0)
    m, t = params
    sampleRate = 20_000  # Hz
    tauSec = (1 / t) / sampleRate

    squaredDiffs = np.square(y - fit_func(x, m, t))
    squaredDiffsFromMean = np.square(y - np.mean(y))
    rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
    print(f"R² = {rSquared}")

    # plot the results
    plt.figure()
    plt.plot(x, y, '.', label="data")
    plt.plot(x, fit_func(x, m, t), '--', label="fitted")
    plt.title("Fitted Power Curve")
    plt.xlabel("r [$\mu m$]")
    plt.ylabel("dV [mV]")
    # plt.show()
    plt.savefig('data/axon_bisc_dist_stim' + f'power_fit_{m} * x^({t}).png')
    # inspect the parameters
    print(f"Y = {m} * x^({t})")
    print(f"Tau = {tauSec * 1e6} µs")


def fit_exponential(x, y, monoExp):
    p0 = (2000, .1, 50)  # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(monoExp, x, y, p0)
    m, t, b = params
    sampleRate = 20_000  # Hz
    tauSec = (1 / t) / sampleRate

    squaredDiffs = np.square(y - monoExp(x, m, t, b))
    squaredDiffsFromMean = np.square(y - np.mean(y))
    rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
    print(f"R² = {rSquared}")

    # plot the results
    plt.figure()
    plt.plot(x, y, '.', label="data")
    plt.plot(x, monoExp(x, m, t, b), '--', label="fitted")
    plt.title("Fitted Exponential Curve")
    plt.show()
    # inspect the parameters
    print(f"Y = {m} * e^(-{t} * x) + {b}")
    print(f"Tau = {tauSec * 1e6} µs")


def fit_linear(x, y, fit_func):
    p0 = (10, 1.3)  # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(fit_func, x, y, p0)
    m, t = params
    sampleRate = 20_000  # Hz
    tauSec = (1 / t) / sampleRate

    squaredDiffs = np.square(y - fit_func(x, m, t))
    squaredDiffsFromMean = np.square(y - np.mean(y))
    rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
    print(f"R² = {rSquared}")

    # plot the results
    plt.figure()
    plt.plot(x, y, '.', label="data")
    plt.plot(x, fit_func(x, m, t), '--', label="fitted")
    plt.title("Fitted linear Curve")
    plt.xlabel("log(r) [$\mu m$]")
    plt.ylabel("log(dV) [mV]")
    # plt.show()
    plt.savefig('data/axon_bisc_dist_stim/' +
                f'linear_fit_{t} * x + {m}_long.png')
    # inspect the parameters
    print(f"Y =  {t} * x + {m}")
    print(f"Tau = {tauSec * 1e6} µs")
