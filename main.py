import sympy as sp


def func(x):
    return 1.62 * (x ** 3) - 8.15 * (x ** 2) + 4.39 * x + 4.29


def graph():
    x = sp.Symbol('x')
    f = func(x)
    f = sp.expand(f)
    sp.plot(f, (x, -2, 5))
    return f


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
        next_x = current_x - (b - current_x)*(f.subs(x, current_x))/(f.subs(x, b) - f.subs(x, current_x))
        while abs(next_x - current_x) > 0.0001:
            current_x = next_x
            next_x = current_x - (b - current_x)*(f.subs(x, current_x))/(f.subs(x, b) - f.subs(x, current_x))
            cnt += 1
    else:
        current_x = b
        next_x = current_x - (current_x - a)*f.subs(x, current_x)/(f.subs(x, current_x) - f.subs(x, a))
        while abs(next_x - current_x):
            current_x = next_x
            next_x = current_x - (current_x - a) * f.subs(x, current_x) / (f.subs(x, current_x) - f.subs(x, a))
            cnt += 1
    print(f'Second method required {cnt} steps')

    # Newton method

    xp = (a+b)/2
    f_derivative = f.diff(x)
    xn = xp - f.subs(x, xp)/f_derivative.subs(x, xp)
    cnt = 0
    while abs(xn - xp) > 0.0001:
        xp = xn
        xn = xp - (f.subs(x, xp)/f_derivative.subs(x, xp))
        cnt += 1

    print(f'Third method required {cnt} steps')
    print(f'{c} {next_x} {xp}\n')


if __name__ == '__main__':
    graph()
    root1, root2, root3 = first_task(-0.6, -0.4), first_task(1, 2), first_task(4, 5)
