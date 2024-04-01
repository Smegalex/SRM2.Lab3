from sympy import pprint, Interval
from sympy.plotting import plot
from sympy.abc import x


def basic_checks(X: list, Y: list):
    if len(Y) != len(X):
        raise (
            ValueError, "Кількість аргументів X не відповідає кількості наданих значень Y")


def multiply_polynomials(polynomial1: dict, polynomial2: dict):
    resulting_polynomial = {}
    for key1, value1 in polynomial1.items():
        key1 = int(key1[2:])  # removing x^ part
        for key2, value2 in polynomial2.items():
            key2 = int(key2[2:])  # removing x^ part
            new_key = f"x^{str(key1 + key2)}"
            new_value = value1 * value2
            if not new_key in resulting_polynomial:
                resulting_polynomial[new_key] = new_value
            else:
                resulting_polynomial[new_key] += new_value
    return resulting_polynomial


def add_polynomials(polynomial1: dict, polynomial2: dict):
    resulting_polynomial = {}
    for polynomial in [polynomial1, polynomial2]:
        for key, value in polynomial.items():
            if not key in resulting_polynomial:
                resulting_polynomial[key] = value
                continue
            resulting_polynomial[key] += value

    return resulting_polynomial


def find_heights(Xs: list) -> dict:
    heights = {}
    for i in range(1, len(Xs)):
        heights[i] = Xs[i]-Xs[i-1]
    return heights


def find_driving_coeficients(driving_coeficients: dict, iteration: int, a: float, b: float, c: float, d: float) -> dict:
    current_A = 0
    current_B = 0
    if a == 0 and c != 0:
        current_A = (-c)/(b)
        current_B = (d)/(b)
    elif a != 0 and c == 0:
        current_A = 0
        current_B = (d-(a*driving_coeficients["B"][iteration-1])) / \
            (b+(a*driving_coeficients["A"][iteration-1]))
    else:
        current_A = (-c)/(b+(a*driving_coeficients["A"][iteration-1]))
        current_B = (d-(a*driving_coeficients["B"][iteration-1])) / \
            (b+(a*driving_coeficients["A"][iteration-1]))

    driving_coeficients["A"][iteration] = current_A
    driving_coeficients["B"][iteration] = current_B
    return driving_coeficients


def solve_der_SLAR(Ys: list, Hs: dict, Qs: dict) -> dict:
    driving_coefficients = {"A": {}, "B": {}}
    for i in range(1, len(Ys)-1):
        a = Hs[i]/6
        b = (Hs[i] + Hs[i+1])/3
        c = Hs[i+1]/6
        if Qs.get(i-1) == 0:
            a = 0
        if Qs.get(i+1) == 0:
            c = 0

        d = ((Ys[i+1]-Ys[i])/Hs[i+1])-((Ys[i]-Ys[i-1])/Hs[i])
        print(d)
        driving_coefficients = find_driving_coeficients(
            driving_coefficients, i, a, b, c, d)
    print(driving_coefficients)
    for j in reversed(range(1, len(Ys)-1)):
        Qs[j] = driving_coefficients["A"][j] * \
            Qs[j+1]+driving_coefficients["B"][j]
    print(Qs)
    return Qs


def calculate_spline(iteration: int, Xs: list, Ys: list, Hs: dict, Qs: dict):
    ind1 = round(Qs[iteration-1]/(6*Hs[iteration]), 5)
    ind2 = round(Qs[iteration]/(6*Hs[iteration]), 5)
    ind3 = round((Ys[iteration-1] / Hs[iteration] -
                 Qs[iteration-1]*Hs[iteration]/6), 5)
    ind4 = round((Ys[iteration]/Hs[iteration] -
                 Qs[iteration] * Hs[iteration]/6), 5)
    spline = ind1*(Xs[iteration]-x)**3+ind2*(x-Xs[iteration-1]
                                             ) ** 3+ind3*(Xs[iteration]-x) + ind4 * (x-Xs[iteration-1])
    print(spline)
    return spline


def find_cubic_spline(Xs: list, Ys: list) -> dict:
    basic_checks(Xs, Ys)
    Qs = {0: 0, 4: 0}
    Hs = find_heights(Xs)

    Qs = solve_der_SLAR(Ys, Hs, Qs)

    splines = {}

    for i in range(1, len(Xs)):
        splines[i] = calculate_spline(i, Xs, Ys, Hs, Qs)
    return splines


def print_n_plot_cubic_spline(splines: dict, Xs: list) -> None:
    SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    ploty = None
    for i in list(splines.keys()):
        print(
            f"S{str(i).translate(SUB)}(x) = ")
        pprint(splines[i])
        print(f"x=[{Xs[i-1]},{Xs[i]}]\n\n")
        if ploty == None:
            ploty = plot(splines[i], (x,
                                      Xs[i-1], Xs[i]), label=f"S{str(i).translate(SUB)}(x)", show=False, legend = True)
            continue

        ploty.extend(plot(splines[i], (x, Xs[i-1], Xs[i]),
                     label=f"S{str(i).translate(SUB)}(x)", show=False))

    return ploty


def find_spline_value(splines: dict, Xs: list, Xstar: float) -> float:
    for i in list(splines.keys()):
        if Interval(Xs[i-1], Xs[i]).contains(Xstar):
            return splines[i].subs(x, Xstar)


def full_cubic_spline_cycle(Xs: list, Ys: list, Xstar: float) -> None:
    splines = find_cubic_spline(Xs, Ys)
    ploty = print_n_plot_cubic_spline(splines, Xs)
    Xstar_value = find_spline_value(splines, Xs, Xstar)
    print(
        f"Кубічний сплайн в точці X* = {Xstar} дорівнює {Xstar_value}.\n\n\n")

    ploty.show()


if __name__ == "__main__":
    Xs = [-3.0, -1.0, 1.0, 3.0, 5.0]
    Ys = [2.81198, 2.3562, 0.7854, 0.32175, 0.19740]
    Xstar = -0.5

    full_cubic_spline_cycle(Xs, Ys, Xstar)
