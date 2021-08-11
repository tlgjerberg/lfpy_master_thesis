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


class Fitting:

    def __init__(self, x, y):

        self.x = x
        self.y = y

    def monoExp(self, x, m, t):
        return m * np.exp(-t * x)

    def powerlaw(self, x, m, t):
        return m * np.power(x, t)

    def linlaw(self, x, m, t):
        return t * x + m

    def curve_fit(self, fit_func="linlaw"):

        fit_functions = {'linlaw': self.linlaw,
                         'monoExp': self.monoExp, 'powerlaw': self.powerlaw}

        self.fit_func = fit_functions[fit_func]

        params, cv = scipy.optimize.curve_fit(
            self.fit_func, x, y, p0)

        self.m, self.t = self.params

    def fit_measure(self):
        squaredDiffs = np.square(y - self.fit_func(x, m, t))
        squaredDiffsFromMean = np.square(y - np.mean(y))
        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
        print(f"R² = {rSquared}")

    def print_fit_func(self):
        if self.fit_func == linlaw:
            print(f"Y =  {t} * x + {m}")
        elif self.fit_func == powerlaw:
            print(f"Y = {m} * x^({t})")
        elif self.fit_func == monoExp:
            print(f"Y = {m} * e^(-{t} * x)")

    def plot_curve(self):

        if self.fit_func == linlaw:
            func_type = 'Linear'
        elif self.fit_func == powerlaw:
            func_type == 'Power'
        elif self.fit_func == monoExp:
            func_type == 'Exponential'

        plt.figure()
        plt.plot(x, y, '.', label="data")
        plt.plot(x, self.fit_func(x, m, t), '--', label="fitted")
        plt.title(f"Fitted {func_type} Curve")
        plt.xlabel("r [$\mu m$]")
        plt.ylabel("dV [mV]")


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
