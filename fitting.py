import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from os.path import join
from sklearn.metrics import r2_score

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

    def print_fit_func(self):
        if self.fit_func == self.linlaw:
            print(f"Y =  {self.t} * x + {self.m}")
        elif self.fit_func == self.powerlaw:
            print(f"Y = {self.m} * x^({self.t})")
        elif self.fit_func == self.monoExp:
            print(f"Y = {self.m} * e^(-{self.t} * x)")

    def fit_measure(self):
        f = self.fit_func(self.x, self.m, self.t)
        squaredDiffs = np.square(self.y - f)
        squaredDiffsFromMean = np.square(self.y - np.mean(self.y))
        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)

        rSquaredAlt = r2_score(self.y, f)

        self.print_fit_func()
        print(f"R² = {rSquared}")
        print(f'R_alt^2 = {rSquaredAlt}')

    def plot_curve(self, xlabel="r [$\mu m$]", ylabel="dV [mV]"):

        if self.fit_func == self.linlaw:
            func_type = 'Linear'
            func_out = f"Y={self.t:.2f}*x+{self.m:.2f}"
        elif self.fit_func == self.powerlaw:
            func_type = 'Power'
            func_out = f"Y={self.m:.2f}*x^({self.t:.2f})"
        elif self.fit_func == self.monoExp:
            func_type = 'Exponential'
            func_out = f"Y={self.m:.2f}*e^(-{self.t:.2f}*x)"

        plt.figure()
        plt.plot(self.x, self.y, '.', label="data")
        plt.plot(self.x, self.fit_func(
            self.x, self.m, self.t), '--', label="fitted")
        plt.title(f"Fitted {func_type} Curve {func_out}")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.savefig(
            join(self.save_folder, f'{func_type}_fit_{func_out}.png'), dpi=300)
