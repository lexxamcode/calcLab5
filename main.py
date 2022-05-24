import sympy as sp
import numpy as np


def func(x):
    return 1.62 * (x ** 3) - 8.15 * (x ** 2) + 4.39 * x + 4.29


def first_graph():
    x = sp.Symbol('x')
    f = func(x)
    f = sp.expand(f)
    sp.plot(f, (x, -2, 5))
    return f


def second_graph():
    x = sp.symbols('x')
    y = sp.symbols('y')
    p1 = sp.plot(sp.sin(x + 1) - 1.2, (x, -2 * sp.pi, 2 * sp.pi), show=False)
    p2 = sp.plot_implicit(2 * x + sp.cos(y) - 2, show=False, line_color='red')
    p1.append(p2[0])
    p1.show()


def first_task(a_arg: float, b_arg: float):
    x = sp.Symbol('x')
    f = func(x)
    f = sp.expand(f)

    a = a_arg
    b = b_arg

    c = (a + b) / 2
    cnt = 0
    while abs(f.subs(x, c)) > 0.0001:
        func_a = f.subs(x, a)
        func_c = f.subs(x, c)
        if (func_a < 0 and func_c > 0) or (func_a > 0 and func_c < 0):
            temp = c
            c = (a + c) / 2
            b = temp
        else:
            temp = c
            c = (b + c) / 2
            a = temp
        cnt += 1
    print(f'First method required {cnt} steps')

    # Horde method

    cnt = 0
    f_derivative = sp.lambdify(x, f)
    if f_derivative(a) > 0:
        current_x = a
        next_x = current_x - (b - current_x) * (f.subs(x, current_x)) / (f.subs(x, b) - f.subs(x, current_x))
        while abs(next_x - current_x) > 0.0001:
            current_x = next_x
            next_x = current_x - (b - current_x) * (f.subs(x, current_x)) / (f.subs(x, b) - f.subs(x, current_x))
            cnt += 1
    else:
        current_x = b
        next_x = current_x - (current_x - a) * f.subs(x, current_x) / (f.subs(x, current_x) - f.subs(x, a))
        while abs(next_x - current_x):
            current_x = next_x
            next_x = current_x - (current_x - a) * f.subs(x, current_x) / (f.subs(x, current_x) - f.subs(x, a))
            cnt += 1
    print(f'Second method required {cnt} steps')

    # Newton method

    xp = (a + b) / 2
    f_derivative = f.diff(x)
    xn = xp - f.subs(x, xp) / f_derivative.subs(x, xp)
    cnt = 0
    while abs(xn - xp) > 0.0001:
        xp = xn
        xn = xp - (f.subs(x, xp) / f_derivative.subs(x, xp))
        cnt += 1

    print(f'Third method required {cnt} steps')
    print(f'{c} {next_x} {xp}\n')


def second_task(x0: float, y0: float):
    x = sp.symbols('x')
    y = sp.symbols('y')

    f1 = sp.sin(x + 1) - y - 1.2
    f2 = 2 * x + sp.cos(y) - 2

    f1_dx = sp.diff(f1, x)
    f1_dy = sp.diff(f1, y)
    f2_dx = sp.diff(f2, x)
    f2_dy = sp.diff(f2, y)

    jacobi_matrix = np.array([[sp.diff(f1, x), sp.diff(f1, y)],
                              [sp.diff(f2, x), sp.diff(f2, y)]])

    # print('\nJacobi matrix:')
    # for column in jacobi_matrix:
    #     print(column)
    # print('\n')

    xn = x0
    yn = y0
    x0 = 0
    count = 0
    eps = 0.0001
    while abs(xn-x0) > eps:
        left_system_matrix = np.array([[f1_dx.subs([(x, xn), (y, yn)]), f1_dy.subs([(x, xn), (y, yn)])],
                                       [f2_dx.subs([(x, xn), (y, yn)]), f2_dy.subs([(x, xn), (y, yn)])]]).astype(float)
        b = np.array([-1*float(f1.subs([(x, xn), (y, yn)])), -1*float(f2.subs([(x, xn), (y, yn)]))])

        # print(f'{left_system_matrix} на {count} итерации')

        temp_matrix = np.copy(left_system_matrix)  # delta_x
        temp_matrix[:, 0] = b
        temp_matrix = temp_matrix.astype(float)
        # print(temp_matrix)
        delta1 = float(np.linalg.det(temp_matrix))
        # print(delta1)

        temp_matrix = np.copy(left_system_matrix)  # delta_y
        temp_matrix[:, 1] = b
        temp_matrix = temp_matrix.astype(float)
        # print(temp_matrix)
        delta2 = float(np.linalg.det(temp_matrix))
        # print(delta2)

        left_system_matrix = left_system_matrix.astype(float)  # delta
        # print(left_system_matrix)
        delta = float(np.linalg.det(left_system_matrix))
        # print(delta)

        delta_x = delta1/delta
        delta_y = delta2/delta

        x0 = xn
        xn = x0 + delta_x
        yn = y0 + delta_y
        count += 1

    print(f'{xn}: {yn} => {count} шагов')


if __name__ == '__main__':
    first_graph()
    root1, root2, root3 = first_task(-0.6, -0.4), first_task(1, 2), first_task(3, 5)
    second_graph()
    second_task(0.5, -0.2)
