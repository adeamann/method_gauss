from numpy import array, zeros, eye, transpose
import sys


def equation_elimination(matrix_a, matrix_b):
    """ Функция прямого хода, применяемая для решения СЛАУ """
    for k in range(n - 1):
        for i in range(k + 1, n):
            if matrix_a[i][k] == 0:
                continue
            factor = matrix_a[i][k] / matrix_a[k][k]
            for j in range(k, n):
                matrix_a[i][j] = matrix_a[i][j] - matrix_a[k][j] * factor
            matrix_b[i] = matrix_b[i] - matrix_b[k] * factor
    return matrix_a, matrix_b


def determinant_elimination(matrix_a):
    """ Функция прямого хода, применяемая для вычисления определителя """
    for k in range(n - 1):
        for i in range(k + 1, n):
            if matrix_a[i][k] == 0:
                continue
            factor = matrix_a[i][k] / matrix_a[k][k]
            for j in range(k, n):
                matrix_a[i][j] = matrix_a[i][j] - matrix_a[k][j] * factor
    return matrix_a


def back_substitution(matrix_a, matrix_b, matrix_x):
    """ Функция обратного хода """
    matrix_x[n - 1] = matrix_b[n - 1] / matrix_a[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        sum_ax = 0
        for j in range(i + 1, n):
            sum_ax += matrix_a[i][j] * matrix_x[j]
        matrix_x[i] = (matrix_b[i] - sum_ax) / matrix_a[i][i]
    return matrix_x


def residual(matrix_a, matrix_b, matrix_x):
    """ Функция вычисления невязки """
    lst = []
    for i in range(n):
        err = 0
        for j in range(n):
            err += matrix_a[i][j] * matrix_x[j]
        err -= matrix_b[i]
        lst.append(err)
    return array(lst)


def det_gauss(matrix_a):
    """ Функция вычисления определителя """
    determination = 1
    for i in range(n):
        determination *= matrix_a[i][i]
    return determination


def reverse_matrix(matrix_a, matrix_e):
    """ Функция нахождения обратной матрицы """
    out = []
    eps = []
    rev_x = zeros(n, float)
    for i in range(n):
        e_col = []
        for j in range(n):
            e_col.append(matrix_e[j][i])
        rev_a, rev_b = equation_elimination(matrix_a, e_col)
        rev_x = back_substitution(rev_a, rev_b, rev_x)
        eps.append(residual(matrix_a, e_col, rev_x))
        out.append(rev_x)
    return transpose(out), array(eps)


if __name__ == '__main__':
    while True:

        a = array([[1.1161, 0.1254, 0.1397, 0.1490],
                   [0.1582, 1.1675, 0.1768, 0.1871],
                   [0.1968, 0.2071, 1.2168, 0.2271],
                   [0.2368, 0.2471, 0.2568, 1.2671]], float)
        n = len(a)

        action = input("Выберите действие: (S) СЛАУ, (D) определитель, (R) обратная матрица, (X) выход: ").lower()

        if action == 's':
            print(f"Порядок матрицы: {n}\n")
            print("Матрица A:\n", a)
            b = array([1.5471, 1.6471, 1.7471, 1.8471], float)
            print("\nВектор свободных членов B:\n", b)
            x = zeros(n, float)
            a_elimination, b_elimination = equation_elimination(a, b)
            print("\nТреугольная матрица:\n", a_elimination)
            print("Преобразованный вектор свободных членов:\n", b_elimination)
            solved_x = back_substitution(a_elimination, b_elimination, x)
            print(f"\nРешение системы: {solved_x}")
            eps_s = residual(a, b, solved_x)
            print(f"Невязка: {eps_s}")

        elif action == 'd':
            print(f"Порядок матрицы: {n}\n")
            print("Матрица A:\n", a)
            a_elimination = determinant_elimination(a)
            print("\nТреугольная матрица:\n", a_elimination)
            d = det_gauss(a_elimination)
            print(f"Определитель матрицы: {d}")

        elif action == 'r':
            print(f"Порядок матрицы: {n}\n")
            print("Матрица A:\n", a)
            e = eye(n)
            print("Единичная матрица E:\n", e)
            rev, eps_r = reverse_matrix(a, e)
            print("Обратная матрица:\n", rev)
            print("Невязка:\n", eps_r)

        elif action == 'x':
            sys.exit()

        else:
            print("Неверный ввод! Пожалуйста, повторите Ваш выбор")
            continue
