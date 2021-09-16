import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from os.path import join

font_params = {
    'font.size': 10,
    'axes.labelsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
}

plt.rcParams.update(**font_params)


class Fitting:

    def __init__(self, save_folder):

        self.save_folder = save_folder

    def monoExp(self, x, m, t):
        return m * np.exp(-t * x)

    def powerlaw(self, x, m, t):
        return m * np.power(x, t)

    def linlaw(self, x, m, t):
        return t * x + m

    def curve_fit(self, x, y, fit_name="linlaw"):

        self.x, self.y = x, y

        fit_functions = {'linlaw': self.linlaw,
                         'monoExp': self.monoExp, 'powerlaw': self.powerlaw}

        self.fit_func = fit_functions[fit_name]

        params, cv = scipy.optimize.curve_fit(self.fit_func, self.x, self.y)

        self.m, self.t = params

    def fit_measure(self):
        squaredDiffs = np.square(y - self.fit_func(x, m, t))
        squaredDiffsFromMean = np.square(y - np.mean(y))
        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
        print(f"R² = {rSquared}")

    def print_fit_func(self):
        if self.fit_func == self.linlaw:
            print(f"Y =  {t} * x + {m}")
        elif self.fit_func == self.powerlaw:
            print(f"Y = {m} * x^({t})")
        elif self.fit_func == self.monoExp:
            print(f"Y = {m} * e^(-{t} * x)")

    def plot_curve(self):

        if self.fit_func == self.linlaw:
            func_type = 'Linear'
        elif self.fit_func == self.powerlaw:
            func_type = 'Power'
        elif self.fit_func == self.monoExp:
            func_type = 'Exponential'

        plt.figure()
        plt.plot(self.x, self.y, '.', label="data")
        plt.plot(self.x, self.fit_func(
            self.x, self.m, self.t), '--', label="fitted")
        plt.title(f"Fitted {func_type} Curve")
        plt.xlabel("r [$\mu m$]")
        plt.ylabel("dV [mV]")

        plt.savefig(join(self.save_folder, func_type + 'fit'))
