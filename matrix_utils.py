# coding=utf8
import extrapolation
from math import pi, sin, e, pow

# lamb = 1.38  # коэф. теплопроводности TODO поменять его, чтобы зависел от т-ры - сделано
u_oc = 293.0  # т-ра для краевых ур-й
alpha = 1.0  # коэф. дл краевых ур-й
f_t = 1.0  # т-ра для краевых условий
k_p = 0.0062  # коэф. поглощения TODO поменять его, чтобы зависел от т-ры
f_0 = 20.0


# параметры - число внутренних точек по осям, шаг сетки (для правой части)
def gen_matr(n_z, n_x, step, lambdas, prev_t):
    n = n_z * n_x
    matrix = [[0.0 for z in xrange(n + 1)] for x in xrange(n)]
    num = 0  # номер уравнения для соответствующей точки
    for j in xrange(1, n_x + 1):
        for i in xrange(1, n_z + 1):
            # выше
            up_i = i
            up_j = j + 1
            if is_inside(up_i, up_j, n_z, n_x):
                pos = number(up_i, up_j, n_z)
                matrix[num][pos] = lam_j_plus(lambdas, i, j)
                # matrix[num][pos] = lam_j_minus_2(i, j, prev_t, n_z)

            # ниже
            down_i = i
            down_j = j - 1
            if is_inside(down_i, down_j, n_z, n_x):
                pos = number(down_i, down_j, n_z)
                matrix[num][pos] = lam_j_minus(lambdas, i, j)
                # matrix[num][pos] = lam_j_plus_2(i, j, prev_t, n_z)

            # слева
            left_i = i - 1
            left_j = j
            if is_inside(left_i, left_j, n_z, n_x):
                pos = number(left_i, left_j, n_z)
                matrix[num][pos] = lam_i_minus(lambdas, i, j)
                # matrix[num][pos] = lam_i_minus_2(i, j, prev_t, n_z)

            # справа
            right_i = i + 1
            right_j = j
            if is_inside(right_i, right_j, n_z, n_x):
                pos = number(right_i, right_j, n_z)
                matrix[num][pos] = lam_i_plus(lambdas, i, j)
                # matrix[num][pos] = lam_i_plus_2(i, j, prev_t, n_z)

            matrix[num][num] = -(lam_j_plus(lambdas, i, j) + lam_j_minus(lambdas, i, j)
                                 + lam_i_plus(lambdas, i, j) + lam_i_minus(lambdas, i, j))
            # matrix[num][num] = -(lam_j_plus_2(i, j, prev_t, n_z) + lam_j_minus_2(i, j, prev_t, n_z)
            #                      + lam_i_plus_2(i, j, prev_t, n_z) + lam_i_minus_2(i, j, prev_t, n_z))

            # правая часть
            x = step * j
            z = step * i
            # matrix[num][n] = sin(pi * x) * sin(pi * z) * step * step
            matrix[num][n] = -f_0 * pow(e, -extrapolation.next_k_p(prev_t[number(i, j, n_z)]) * z) * step * step
            # matrix[num][n] = -f_0 * pow(e, -k_p * z) * step * step

            # добавление краевых условий
            if is_on_left_bound(i):
                matrix[num][num] += lam_i_minus(lambdas, i, j)
                matrix[num][n] -= lam_i_minus(lambdas, i, j) * step * f_t / lambdas[0][j]
            if is_on_right_bound(i, n_z):
                matrix[num][num] += lam_i_plus(lambdas, i, j) * lambdas[n_z + 1][j] / (lambdas[n_z + 1][j] + alpha * step)
                matrix[num][n] -= lam_i_plus(lambdas, i, j) * alpha * step * u_oc / (lambdas[n_z + 1][j] + alpha * step)
            if is_on_lower_bound(j):
                matrix[num][num] += lam_j_minus(lambdas, i, j) * lambdas[i][0] / (lambdas[i][0] - alpha * step)
                matrix[num][n] += lam_j_minus(lambdas, i, j) * alpha * step * u_oc / (lambdas[i][0] - alpha * step)
            if is_on_upper_bound(j, n_x):
                matrix[num][num] += lam_j_plus(lambdas, i, j) * lambdas[i][n_x + 1] / (lambdas[i][n_x + 1] + alpha * step)
                matrix[num][n] -= lam_j_plus(lambdas, i, j) * alpha * step * u_oc / (lambdas[i][n_x + 1] + alpha * step)

            num += 1

    # print output_matrix(matrix)
    return matrix


# проверка точки на принадлежность к внутренней
def is_inside(i, j, n_z, n_x):
    if 1 <= i <= n_z and 1 <= j <= n_x:
        return True
    else:
        return False


# номер в столбце неизвестных для соответствующей точки
def number(i, j, n_z):
    return (j - 1) * n_z + i - 1


# проверка на принадлежность к нижней границе
def is_on_lower_bound(j):
    return j == 1


# проверка на принадлежность к верхней границе
def is_on_upper_bound(j, n_x):
    return j == n_x


# проверка на принадлежность к левой границе
def is_on_left_bound(i):
    return i == 1


# проверка на принадлежность к правой границе
def is_on_right_bound(i, n_z):
    return i == n_z


# пересчитывает матрицу коэффициентов теплопроводности
def recalc_lambda_coef(coefs, temps):
    m = len(coefs)  # перебираем по z
    n = len(coefs[0])  # перебираем по x
    n_z = m - 2
    for i in xrange(m):
        for j in xrange(n):
            if i != 0 and i != m - 1 and j != 0 and j != n - 1:
                coefs[i][j] = lambda_coef(temps[number(i, j, n_z)])
            else:
                coefs[i][j] = lambda_coef(20.0)
    return coefs


# формула расчёта коэффициента теплопроводности от т-ры
def lambda_coef(t):
    return 0.134 * 0.1 * (1 + 4.35 * 0.0001 * t)


# полуцелый коэф. теплопроводности
def lam_j_minus(lambdas, i, j):
    return (lambdas[i][j - 1] + lambdas[i][j]) / 2


def lam_j_plus(lambdas, i, j):
    return (lambdas[i][j + 1] + lambdas[i][j]) / 2


def lam_i_minus(lambdas, i, j):
    return (lambdas[i - 1][j] + lambdas[i][j]) / 2


def lam_i_plus(lambdas, i, j):
    return (lambdas[i + 1][j] + lambdas[i][j]) / 2


def lam_j_minus_2(i, j, t, n_z):
    p = number(i, j - 1, n_z)
    if not (0 <= p < len(t)):
        p = 20.0
    else:
        p = t[number(i, j - 1, n_z)]
    return lambda_coef((p + t[number(i, j, n_z)]) / 2)


def lam_j_plus_2(i, j, t, n_z):
    p = number(i, j + 1, n_z)
    if not (0 <= p < len(t)):
        p = 20.0
    else:
        p = t[number(i, j + 1, n_z)]
    return lambda_coef((p + t[number(i, j, n_z)]) / 2)


def lam_i_minus_2(i, j, t, n_z):
    p = number(i - 1, j, n_z)
    if not (0 <= p < len(t)):
        p = 20.0
    else:
        p = t[number(i - 1, j, n_z)]
    return lambda_coef((p + t[number(i, j, n_z)]) / 2)


def lam_i_plus_2(i, j, t, n_z):
    p = number(i + 1, j, n_z)
    if not (0 <= p < len(t)):
        p = 20.0
    else:
        p = t[number(i + 1, j, n_z)]
    return lambda_coef((p + t[number(i, j, n_z)]) / 2)


# генерации массивов коэффициентов диагональных элементов для метода прогонки (промежуточный слой)
# i должен быть внутри прямоугольника
def gen_pr_matr_prom(i, n_z, n_x, step, t_step, lambdas, y):
    a = [0.0] * n_x  # под главной диагональю
    b = [0.0] * n_x  # главная диагональ
    c = [0.0] * (n_x - 1)  # над главной диагональю
    f = [0.0] * n_x  # столбец правых частей
    z = i * step  # координата точки в прямоугольнике
    for j in xrange(2, n_x + 1):
        a[j - 1] = lam_j_minus(lambdas, i, j)
    for j in xrange(1, n_x + 1):
        b[j - 1] = -(step * step / t_step + lam_j_plus(lambdas, i, j) + lam_j_minus(lambdas, i, j))
        if j == 1:  # кр. условие снизу
            b[j - 1] += lam_j_minus(lambdas, i, j) * lambdas[i][0] / (lambdas[i][0] - alpha * step)
        if j == n_x:  # кр. условие сверху
            b[j - 1] += lam_j_plus(lambdas, i, j) * lambdas[i][n_x + 1] / (lambdas[i][n_x + 1] + alpha * step)
    for j in xrange(1, n_x):
        c[j - 1] = lam_j_plus(lambdas, i, j)
    for j in xrange(1, n_x + 1):
        # f[j - 1] = (-f_0 * pow(e, -k_p * z) / 2.0 - y[number(i, j, n_z)] / t_step) * step * step
        f[j - 1] = (-f_0 * pow(e, -extrapolation.next_k_p(y[number(i, j, n_z)]) * z) / 2.0
                    - y[number(i, j, n_z)] / t_step) * step * step
        if j == 1:  # кр. условие снизу
            f[j - 1] += lam_j_minus(lambdas, i, j) * alpha * step * u_oc / (lambdas[i][0] - alpha * step)
        if j == n_x:  # кр. условие сверху
            f[j - 1] -= lam_j_plus(lambdas, i, j) * alpha * step * u_oc / (lambdas[i][n_x + 1] + alpha * step)

    return a, b, c, f


# генерации массивов коэффициентов диагональных элементов для метода прогонки (новый слой)
# j должен быть внутри прямоугольника
def gen_pr_matr_new(j, n_z, n_x, step, t_step, lambdas, y):
    a = [0.0] * n_z  # под главной диагональю
    b = [0.0] * n_z  # главная диагональ
    c = [0.0] * (n_z - 1)  # над главной диагональю
    f = [0.0] * n_z  # столбец правых частей
    for i in xrange(2, n_z + 1):
        a[i - 1] = lam_i_minus(lambdas, i, j)
    for i in xrange(1, n_z + 1):
        b[i - 1] = -(step * step / t_step + lam_i_plus(lambdas, i, j) + lam_i_minus(lambdas, i, j))
        if i == 1:  # кр. условие слева
            b[i - 1] += lam_i_minus(lambdas, i, j)
        if i == n_z:  # кр. условие справа
            b[i - 1] += lam_i_plus(lambdas, i, j) * lambdas[n_z + 1][j] / (lambdas[n_z + 1][j] + alpha * step)
    for i in xrange(1, n_z):
        c[i - 1] = lam_i_plus(lambdas, i, j)
    for i in xrange(1, n_z + 1):
        z = i * step  # координата точки в прямоугольнике
        # f[i - 1] = (-f_0 * pow(e, -k_p * z) / 2.0 - y[number(i, j, n_z)] / t_step) * step * step
        f[i - 1] = (-f_0 * pow(e, -extrapolation.next_k_p(y[number(i, j, n_z)]) * z) / 2.0
                    - y[number(i, j, n_z)] / t_step) * step * step
        if i == 1:  # кр. условие слева
            f[i - 1] -= lam_i_minus(lambdas, i, j) * step * f_t / lambdas[0][j]
        if i == n_z:  # кр. условие справа
            f[i - 1] -= lam_i_plus(lambdas, i, j) * alpha * step * u_oc / (lambdas[n_z + 1][j] + alpha * step)

    return a, b, c, f


# вывод
def output_matrix(a):
    for i in xrange(len(a)):
        for j in xrange(len(a[i])):
            print "{:4}".format(a[i][j]),
            # print a[i][j],
        print
