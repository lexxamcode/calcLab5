from matplotlib import pyplot as plt
import numpy as np


def show_polynomial(a: float, b: float,
                    cube_coefficient: float, quad_coefficient: float, x_coefficient: float, coefficient: float):
    plt.figure(figsize=(5, 3), dpi=80)  # Window create
    ax = plt.subplot()

    ax.spines['right'].set_color('none')  # Do not show right and top borders
    ax.spines['top'].set_color('none')

    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data', 0))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data', 0))

    x = np.linspace(a, b, 256, endpoint=True)
    func = (cube_coefficient*(x**3)) + (quad_coefficient*(x**2)) + (x_coefficient*x) + coefficient

    plt.plot(x, func, color='blue')  # Adds the graph to our plot

    plt.legend(loc='upper left', frameon=False)
    plt.grid()
    plt.show()


if __name__ == '__main__':
    show_polynomial(-2, 5, 1.62, -8.15, 4.39, 4.29)
