# coding=utf8
from math import pi, sin, e, pow

lamb = 1.38  # коэф. теплопроводности
u_oc = 293  # т-ра для краевых ур-й
alpha = 1  # коэф. дл краевых ур-й
f_t = 1  # т-ра для краевых условий
k_p = 0.0062  # коэф. поглощения
f_0 = 1


# параметры - число внутренних точек по осям, шаг сетки (для правой части)
def gen_matr(n_z, n_x, step):
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
                matrix[num][pos] = -1.0

            # ниже
            down_i = i
            down_j = j - 1
            if is_inside(down_i, down_j, n_z, n_x):
                pos = number(down_i, down_j, n_z)
                matrix[num][pos] = -1.0

            # слева
            left_i = i - 1
            left_j = j
            if is_inside(left_i, left_j, n_z, n_x):
                pos = number(left_i, left_j, n_z)
                matrix[num][pos] = -1.0

            # справа
            right_i = i + 1
            right_j = j
            if is_inside(right_i, right_j, n_z, n_x):
                pos = number(right_i, right_j, n_z)
                matrix[num][pos] = -1.0

            matrix[num][num] = 4.0

            # правая часть
            x = step * j
            z = step * i
            # matrix[num][n] = sin(pi * x) * sin(pi * z) * step * step
            matrix[num][n] = f_0 * pow(e, -k_p * z) * step * step

            # добавление краевых условий
            if is_on_left_bound(i):
                matrix[num][num] -= 1
                matrix[num][n] += f_t * step / lamb
            if is_on_right_bound(i, n_z):
                matrix[num][num] -= lamb / (lamb + alpha * step)
                matrix[num][n] += alpha * step * u_oc / (lamb + alpha * step)
            if is_on_lower_bound(j):
                matrix[num][num] -= lamb / (lamb - alpha * step)
                matrix[num][n] += alpha * step * u_oc / (lamb - alpha * step)
            if is_on_upper_bound(j, n_x):
                matrix[num][num] -= lamb / (lamb + alpha * step)
                matrix[num][n] += alpha * step * u_oc / (lamb + alpha * step)

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


# вывод
def output_matrix(a):
    for i in xrange(len(a)):
        for j in xrange(len(a[i])):
            print "{:4}".format(a[i][j]),
            # print a[i][j],
        print
