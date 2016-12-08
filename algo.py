# coding=utf8

eps = 0.00001  # точность выбора опорного элемента
eps_conv = 0.001  # точность сходимости


# метод Гаусса
# a - исходная двумерная матрица коэффициентов
def gauss(a):
    # Прямой
    # ход

    n = len(a)  # число строк
    m = len(a[0])  # число столбцов

    row = 0
    col = 0
    while row < n and col < m:
        # находим строку с опорным элементом (опорный элемент - наибольший по модулю)
        sel = row
        for i in xrange(row + 1, n):
            if abs(a[i][col]) > abs(a[sel][col]):
                sel = i
        if abs(a[sel][col]) < eps:
            col += 1
            continue

        # меняем текущую и опорные строки местами
        for j in xrange(col, m):
            temp = a[sel][j]
            a[sel][j] = a[row][j]
            a[row][j] = temp

        # зануляем столбец col для строк, ниже текущей = опорной, вычитая из каждой опорную, умнож. на коэф.
        for i in xrange(row + 1, n):
            c = a[i][col] / a[row][col]
            for j in xrange(col, m):
                a[i][j] -= a[row][j] * c

        col += 1
        row += 1

    # Обратный
    # ход

    ans = [0] * n

    ans[n - 1] = a[n - 1][m - 1] / a[n - 1][m - 2]
    for i in xrange(n - 2, -1, -1):
        s = a[i][m - 1]  # правая часть
        for j in xrange(i + 1, m - 1):
            s -= a[i][j] * ans[j]
        ans[i] = s / a[i][i]

    return a, ans


# метод прогонки
# b - диагональ под главной(n-1), d - диагональ над главной(n-1), c - главная диагональ(n), r - правые части(n)
# b - [1..n-1], d - [0..n-2], c - [0..n-1], r - [0..n-1]
def tdma(b, d, c, r):
    # преобразовываем к дробному виду
    b, c, d, r = map(lambda k_list: map(float, k_list), (b, c, d, r))

    # прогоночные коэффициенты
    delta = [-d[0] / c[0]]
    lamd = [r[0] / c[0]]
    n = len(r)
    x = [0] * n  # искомые неизвестные

    # прямой ход
    for i in xrange(1, n):
        if i != n - 1:
            delta.append(-d[i] / (c[i] + b[i] * delta[i - 1]))
        else:
            delta.append(0)
        lamd.append((r[i] - b[i] * lamd[i - 1]) / (c[i] + b[i] * delta[i - 1]))

    # обратный ход
    x[n - 1] = lamd[n - 1]
    for i in xrange(n - 2, -1, -1):
        x[i] = delta[i] * x[i + 1] + lamd[i]

    return x


# проверка сходимости
def is_conv(y_old, y_new):
    max_diff = -1
    for i in xrange(len(y_old)):
        if abs(y_old[i] - y_new[i]) > max_diff:
            max_diff = abs(y_old[i] - y_new[i])
    if max_diff >= eps_conv:
        # print max_diff
        return False
    return True