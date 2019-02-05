#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#         Реализация двух и трехмерных конечных элементов
###################################################################

import math
from abc import abstractmethod
from numpy.linalg import solve, LinAlgError
from numpy import array
from numpy import zeros
from numpy import identity
from numpy.linalg import det
from numpy.linalg import inv
from numpy.linalg import norm
from numpy import cross
from core.fem_error import TFEMException
from core.fem_defs import eps
from core.fem_params import TFEMParams


# Векторное произведение a и b
def cross_product(a, b):
    v = cross(a, b)
    return v / norm(v)


# Построение вектора для заданных координат
def vector(x0, x1):
    v = array(x0) - array(x1)
    return v / norm(v)


# Матрица преобразования в локальную систему координат
def create_transform_matrix(x):
    v_x = vector(x[1], x[0])
    v_z = cross_product(vector(x[1], x[0]), vector(x[2], x[0]))
    v_y = cross_product(v_z, v_x)
    return array([v_x, v_y, v_z])


# Абстрактный базовый класс, описывающий конечный элемент (КЭ)
class TFE:
    def __init__(self):
        self.size = 0               # Размерность КЭ
        self.params = TFEMParams()  # Параметры (упругие свойства, толщина (для оболочек и пластин) и т.п.)
        self.x = []                 # Координаты вершин КЭ
        self.K = []                 # Локальная матрица жесткости
        self.M = []                 # ... масс
        self.C = []                 # ... демпфирования
        self.load = []              # Локальный столбец нагрузки
        self.a = []                 # Коэффициенты функций форм

    # Задание параметров
    def set_params(self, p):
        self.params = p

    # Задание координат
    def set_coord(self, x):
        self.x = array(x)
        self._create()

    # Вычисление функций форм КЭ
    @abstractmethod
    def _create(self):
        raise NotImplementedError('Method TFE._create() is pure virtual')

    # Матрица упругих свойств
    @abstractmethod
    def _elastic_matrix(self):
        raise NotImplementedError('Method TFE._elastic_matrix() is pure virtual')

    # Вычисления стандартных результатов КЭ
    @abstractmethod
    def calc(self, u):
        raise NotImplementedError('Method TFE.calc() is pure virtual')

    # Вычисление матриц жесткости, масс и демпфирования
    @abstractmethod
    def generate(self, v_load, s_load, is_static=True):
        raise NotImplementedError('Method TFE.generate() is pure virtual')


# Линейный (двухузловой) одномерный КЭ
class TFE1D2(TFE):
    def __init__(self):
        super().__init__()
        self.size = 2

    def _length(self):
        return math.fabs(self.x[1][0] - self.x[0][0])

    def _create(self):
        if self._length() == 0.0:
            raise TFEMException('incorrect_fe_err')
        self.a = zeros((self.size, self.size))
        self.a[0][0] = self.x[1][0]/(self.x[1][0] - self.x[0][0])
        self.a[0][1] = -1.0/(self.x[1][0] - self.x[0][0])
        self.a[1][0] = self.x[0][0]/(self.x[0][0] - self.x[1][0])
        self.a[1][1] = -1.0/(self.x[0][0] - self.x[1][0])

    def calc(self, u):
        res = zeros((self.size, self.size))
        res[0][0] = res[0][1] = u[0] * self.a[0][1] + u[1] * self.a[1][1]
        res[1][0] = res[1][1] = self.params.e[0] * (u[0] * self.a[0][1] + u[1] * self.a[1][1])
        return res

    def _elastic_matrix(self):
        return array([
            [1.0, -1.0],
            [-1.0, 1.0]
        ]) * self.params.e[0]

    def generate(self, v_load, s_load, is_static=True):
        self.K = self._elastic_matrix() / self._length() * self.params.thickness
        self.load = array([1 / 2, 1 / 2]) * self._length * self.params.thickness + array([1, 1]) * self.params.thickness # Уточнить вычисление поверхностной нагрузки
        if not is_static:
            m = array([
                [1.0, -1.0],
                [-1.0, 1.0]
            ])
            self.M = 1.0/6.0 * m * self._length() * self.params.thickness * self.params.density
            self.C = 1.0/6.0 * m * self._length() * self.params.thickness * self.params.damping


# Линейный (трехузловой) треугольный КЭ
class TFE2D3(TFE):
    def __init__(self):
        super().__init__()
        self.size = 3
        self._xi = [0, 1 / 2, 1 / 2]
        self._eta = [1 / 2, 0, 1 / 2]
        self._w = [1 / 6, 1 / 6, 1 / 6]

    def _elastic_matrix(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.params.m[0], 0.0],
            [self.params.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - self.params.m[0])]
        ]) * self.params.e[0]/(1.0 - self.params.m[0]**2)
        return d

    def _square(self):
        a = math.sqrt((self.x[0][0] - self.x[1][0])**2 + (self.x[0][1] - self.x[1][1])**2)
        b = math.sqrt((self.x[0][0] - self.x[2][0])**2 + (self.x[0][1] - self.x[2][1])**2)
        c = math.sqrt((self.x[2][0] - self.x[1][0])**2 + (self.x[2][1] - self.x[1][1])**2)
        p = 0.5 * (a + b + c)
        return math.sqrt(p * (p - a) * (p - b) * (p - c))

    def _create(self):
        det0 = self.x[2][1] * self.x[1][0] - self.x[2][1] * self.x[0][0] - self.x[0][1] * self.x[1][0] - \
               self.x[1][1] * self.x[2][0] + self.x[1][1] * self.x[0][0] + self.x[0][1] * self.x[2][0]
        if math.fabs(det0) < eps:
            raise TFEMException('incorrect_fe_err')
        index = [[2, 1], [0, 2], [1, 0]]
        self.a = zeros((self.size, self.size))
        for i in range(0, self.size):
            det1 = self.x[index[i][0]][1] * self.x[index[i][1]][0] - self.x[index[i][1]][1] * self.x[index[i][0]][0]
            det2 = self.x[index[i][1]][1] - self.x[index[i][0]][1]
            det3 = self.x[index[i][0]][0] - self.x[index[i][1]][0]
            self.a[i][0] = det1/det0
            self.a[i][1] = det2/det0
            self.a[i][2] = det3/det0

    def calc(self, u):
        res = zeros((6, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([3, 2 * self.size])
            for j in range(0, self.size):
                b[0][j * 2 + 0] = b[2][j * 2 + 1] = self.a[j][1]    # dx
                b[1][j * 2 + 1] = b[2][j * 2 + 0] = self.a[j][2]    # dy
            e = b.dot(u)
            s = self._elastic_matrix().dot(e)
            for j in range(0, 3):
                res[j][i] += e[j]
                res[j + 3][i] += s[j]
        return res

    # Изопараметрические функции формы и их производные
    def _shape(self, i):
        return array([1 - self._xi[i] - self._eta[i], self._xi[i], self._eta[i]])

    def _shape_dxi(self, i):
        return array([-1, 1, 0])

    def _shape_deta(self, i):
        return array([-1, 0, 1])

    # Формирование локальной матрицы жесткости
    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((2 * self.size, 2 * self.size))
        self.load = zeros(2 * self.size)
        if not is_static:
            self.M = zeros((2 * self.size, 2 * self.size))
            self.C = zeros((2 * self.size, 2 * self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[0,1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Изопараметрические функции формы и их производные
            # Матрица Якоби
            jacobi = array([
                [sum(self._shape_dxi(i) * self.x[:, 0]), sum(self._shape_dxi(i) * self.x[:, 1])],
                [sum(self._shape_deta(i) * self.x[:, 0]), sum(self._shape_deta(i) * self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inv_jacobi = inv(jacobi)
            shape_dx = inv_jacobi[0, 0] * self._shape_dxi(i) + inv_jacobi[0, 1] * self._shape_deta(i)
            shape_dy = inv_jacobi[1, 0] * self._shape_dxi(i) + inv_jacobi[1, 1] * self._shape_deta(i)
            # Матрицы градиентов и функций форм
            b = zeros([3, 2 * self.size])
            n = zeros([2, 2 * self.size])
            for j in range(0, self.size):
                b[0][2 * j + 0] = b[2][2 * j + 1] = shape_dx[j]
                b[1][2 * j + 1] = b[2][2 * j + 0] = shape_dy[j]
                n[1][2 * j + 1] = n[0][2 * j + 0] = self._shape(i)[j]
            # Вычисление компонент локальной матрицы жесткости
            self.K += (b.conj().transpose().dot(self._elastic_matrix()).dot(b) * self.params.thickness *
                       abs(jacobian) * self._w[i])
            self.load += n.conj().transpose().dot(v_load) * self.params.thickness * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += (self._shape(i).conj().transpose().dot(self._shape(i)) * self.params.thickness *
                           abs(jacobian) * self._w[i] * self.params.density)
                self.C += (self._shape(i).conj().transpose().dot(self._shape(i)) * self.params.thickness *
                           abs(jacobian) * self._w[i] * self.params.damping)


# Квадратичный (шестиузловой) треугольный КЭ
class TFE2D6(TFE2D3):
    def __init__(self):
        super().__init__()
        self.size = 6
        # Параметры квадратур Гаусса
        self._xi = [1 / 3, 0, 1 / 2, 1 / 2, 1, 0, 0]
        self._eta = [1 / 3, 1 / 2, 0, 1 / 2, 0, 1, 0]
        self._w = [27 / 120, 8 / 120, 8 / 120, 8 / 120, 3 / 120, 3 / 120, 3 / 120]

    def _create(self):
        if self._square() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][0] * self.x[i][1]
                a[i][4] = self.x[i][0] ** 2
                a[i][5] = self.x[i][1] ** 2
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    # Изопараметрические функции формы и их производные
    def _shape(self, i):
        s = array([1 - self._xi[i] - self._eta[i], self._xi[i], self._eta[i]])
        return array([s[0] * (2 * s[0] - 1), s[1] * (2 * s[1] - 1), s[2] * (2 * s[2] - 1), 4 * s[0] * s[1],
                      4 * s[1] * s[2], 4 * s[0] * s[2]])

    def _shape_dxi(self, i):
        return array([-3 + 4 * self._xi[i] + 4 * self._eta[i], 4 * self._xi[i] - 1, 0,
                      -8 * self._xi[i] + 4 - 4 * self._eta[i], 4 * self._eta[i], -4 * self._eta[i]])

    def _shape_deta(self, i):
        return array([-3 + 4 * self._xi[i] + 4 * self._eta[i], 0, 4 * self._eta[i] - 1, -4 * self._xi[i],
                      4 * self._xi[i], -8 * self._eta[i] + 4 - 4 * self._xi[i]])

    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((2 * self.size, 2 * self.size))
        self.load = zeros(2 * self.size)
        if not is_static:
            self.M = zeros((2 * self.size, 2 * self.size))
            self.C = zeros((2 * self.size, 2 * self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[0,1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Матрица Якоби
            jacobi = array([
                [sum(self._shape_dxi(i) * self.x[:, 0]), sum(self._shape_dxi(i) * self.x[:, 1])],
                [sum(self._shape_deta(i) * self.x[:, 0]), sum(self._shape_deta(i) * self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inv_jacobi = inv(jacobi)
            shape_dx = inv_jacobi[0, 0] * self._shape_dxi(i) + inv_jacobi[0, 1] * self._shape_deta(i)
            shape_dy = inv_jacobi[1, 0] * self._shape_dxi(i) + inv_jacobi[1, 1] * self._shape_deta(i)
            # Матрицы градиентов и функций форм
            b = zeros([3, 2 * self.size])
            n = zeros([2, 2 * self.size])
            for j in range(0, self.size):
                b[0][2 * j + 0] = b[2][2 * j + 1] = shape_dx[j]
                b[1][2 * j + 1] = b[2][2 * j + 0] = shape_dy[j]
                n[0][2 * j + 0] = n[1][2 * j + 1] = self._shape(i)[j]
            # Вычисление компонент локальной матрицы жесткости
            self.K += (b.conj().transpose().dot(self._elastic_matrix()).dot(b) * self.params.thickness *
                       abs(jacobian) * self._w[i])
            self.load += n.conj().transpose().dot(v_load) * self.params.thickness * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += (self._shape(i).conj().transpose().dot(self._shape(i)) * self.params.thickness *
                           abs(jacobian) * self._w[i] * self.params.density)
                self.C += (self._shape(i).conj().transpose().dot(self._shape(i)) * self.params.thickness *
                           abs(jacobian) * self._w[i] * self.params.damping)

    def calc(self, u):
        res = zeros((6, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([3, 2 * self.size])
            for j in range(0, self.size):
                dx = self.a[j][1] + self.a[j][3] * self.x[i][1] + 2 * self.a[j][4] * self.x[i][0]
                dy = self.a[j][2] + self.a[j][3] * self.x[i][0] + 2 * self.a[j][5] * self.x[i][1]
                b[0][j * 2 + 0] = b[2][j * 2 + 1] = dx
                b[1][j * 2 + 1] = b[2][j * 2 + 0] = dy
            e = b.dot(u)
            s = self._elastic_matrix().dot(e)
            for j in range(0, 3):
                res[j][i] += e[j]
                res[j + 3][i] += s[j]
        return res


# Билинейный четырехузловой двумерный КЭ
class TFE2D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4
        # Параметры квадратур Гаусса
        self._xi = [-0.57735027, -0.57735027, 0.57735027, 0.57735027]
        self._eta = [-0.57735027, 0.57735027, -0.57735027, 0.57735027]
        self._w = [1.0, 1.0, 1.0, 1.0]

    def _square(self):
        return math.sqrt((self.x[0][0] - self.x[1][0])**2 + (self.x[0][1] - self.x[1][1])**2)

    def _create(self):
        if self._square() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][0] * self.x[i][1]
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    def _elastic_matrix(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.params.m[0], 0.0],
            [self.params.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - self.params.m[0])]
        ]) * self.params.e[0]/(1.0 - self.params.m[0]**2)
        return d

    def calc(self, u):
        res = zeros((6, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([3, 2 * self.size])
            for j in range(0, self.size):
                b[0][j * 2 + 0] = b[2][j * 2 + 1] = self.a[j][1] + self.a[j][3] * self.x[i][1]  # dx
                b[1][j * 2 + 1] = b[2][j * 2 + 0] = self.a[j][2] + self.a[j][3] * self.x[i][0]  # dy
            e = b.dot(u)
            s = self._elastic_matrix().dot(e)
            for j in range(0, 3):
                res[j][i] += e[j]
                res[j + 3][i] += s[j]
        return res

    # Изопараметрические функции формы и их производные
    def _shape(self, i):
        return array([
            0.25 * (1.0 - self._xi[i]) * (1.0 - self._eta[i]),
            0.25 * (1.0 + self._xi[i]) * (1.0 - self._eta[i]),
            0.25 * (1.0 + self._xi[i]) * (1.0 + self._eta[i]),
            0.25 * (1.0 - self._xi[i]) * (1.0 + self._eta[i])
        ])

    def _shape_dxi(self, i):
        return array([
            -0.25 * (1.0 - self._eta[i]),
            0.25 * (1.0 - self._eta[i]),
            0.25 * (1.0 + self._eta[i]),
            -0.25 * (1.0 + self._eta[i])
        ])

    def _shape_deta(self, i):
        return array([
            -0.25 * (1.0 - self._xi[i]),
            -0.25 * (1.0 + self._xi[i]),
            0.25 * (1.0 + self._xi[i]),
            0.25 * (1.0 - self._xi[i])
        ])

    # Формирование локальной матрицы жесткости
    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((2 * self.size, 2 * self.size))
        self.load = zeros(2 * self.size)
        if not is_static:
            self.M = zeros((2 * self.size, 2 * self.size))
            self.C = zeros((2 * self.size, 2 * self.size))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Матрица Якоби
            jacobi = array([
                [sum(self._shape_dxi(i) * self.x[:, 0]), sum(self._shape_dxi(i) * self.x[:, 1])],
                [sum(self._shape_deta(i) * self.x[:, 0]), sum(self._shape_deta(i) * self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            if jacobian <= 0:
                print('Jacobian error %f' % jacobian)
            inv_jacobi = inv(jacobi)
            shape_dx = inv_jacobi[0, 0] * self._shape_dxi(i) + inv_jacobi[0, 1] * self._shape_deta(i)
            shape_dy = inv_jacobi[1, 0] * self._shape_dxi(i) + inv_jacobi[1, 1] * self._shape_deta(i)
            # Изопараметрические матрицы градиентов и функций форм
            b = zeros([3, 2 * self.size])
            n = zeros([2, 2 * self.size])
            for j in range(0, self.size):
                b[0][j * 2 + 0] = b[2][j * 2 + 1] = shape_dx[j]
                b[1][j * 2 + 1] = b[2][j * 2 + 0] = shape_dy[j]
                n[0][2 * j + 0] = n[1][2 * j + 1] = self._shape(i)[j]

            self.K += b.conj().transpose().dot(self._elastic_matrix()).dot(b) * abs(jacobian) * self._w[i] *  \
                      self.params.thickness
            self.load += n.conj().transpose().dot(v_load) * self.params.thickness * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += (self._shape(i).conj().transpose().dot(self._shape(i)) * self.params.thickness *
                           abs(jacobian) * self._w[i] * self.params.density)
                self.C += (self._shape(i).conj().transpose().dot(self._shape(i)) * self.params.thickness *
                           abs(jacobian) * self._w[i] * self.params.damping)


# Линейный (четырехузловой) тетраэдральный КЭ
class TFE3D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4
        # Параметры квадратур Гаусса
        self._xi = [1 / 4, 1 / 2, 1 / 6, 1 / 6, 1 / 6]
        self._eta = [1 / 4, 1 / 6, 1 / 2, 1 / 6, 1 / 6]
        self._psi = [1 / 4, 1 / 6, 1 / 6, 1 / 2, 1 / 6]
        self._w = [-4 / 30, 9 / 120, 9 / 120, 9 / 120, 9 / 120]

    def _elastic_matrix(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.params.m[0]/(1.0 - self.params.m[0]), self.params.m[0]/(1.0 - self.params.m[0]), 0.0, 0.0, 0.0],
            [self.params.m[0]/(1.0 - self.params.m[0]), 1.0, self.params.m[0]/(1.0 - self.params.m[0]), 0.0, 0.0, 0.0],
            [self.params.m[0]/(1.0 - self.params.m[0]), self.params.m[0]/(1.0 - self.params.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * self.params.m[0])/(1.0 - self.params.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * self.params.m[0])/(1.0 - self.params.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * self.params.m[0])/(1.0 - self.params.m[0])],
        ]) * self.params.e[0] * (1.0 - self.params.m[0])/(1.0 + self.params.m[0])/(1.0 - 2.0 * self.params.m[0])
        return d

    def _volume(self):
        a = (self.x[1][0] - self.x[0][0]) * (self.x[2][1] - self.x[0][1]) * (self.x[3][2] - self.x[0][2]) + \
            (self.x[3][0] - self.x[0][0]) * (self.x[1][1] - self.x[0][1]) * (self.x[2][2] - self.x[0][2]) + \
            (self.x[2][0] - self.x[0][0]) * (self.x[3][1] - self.x[0][1]) * (self.x[1][2] - self.x[0][2])
        b = (self.x[3][0] - self.x[0][0]) * (self.x[2][1] - self.x[0][1]) * (self.x[1][2] - self.x[0][2]) + \
            (self.x[2][0] - self.x[0][0]) * (self.x[1][1] - self.x[0][1]) * (self.x[3][2] - self.x[0][2]) + \
            (self.x[1][0] - self.x[0][0]) * (self.x[3][1] - self.x[0][1]) * (self.x[2][2] - self.x[0][2])
        return math.fabs(a - b)/6.0

    def _create(self):
        if self._volume() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][2]
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([6, 3 * self.size])
            for j in range(0, self.size):
                b[0][j * 3 + 0] = b[3][j * 3 + 1] = b[5][j * 3 + 2] = self.a[j][1]  # dx
                b[1][j * 3 + 1] = b[3][j * 3 + 0] = b[4][j * 3 + 2] = self.a[j][2]  # dy
                b[2][j * 3 + 2] = b[4][j * 3 + 1] = b[5][j * 3 + 0] = self.a[j][3]  # dz
            e = b.dot(u)
            s = self._elastic_matrix().dot(e)
            for j in range(0, 6):
                res[j][i] += e[j]
                res[j + 6][i] += s[j]
        return res

    # Изопараметрические функции формы и их производные
    def _shape(self, i):
        return array([1 - self._xi[i] - self._eta[i] - self._psi[i], self._xi[i], self._eta[i], self._psi[i]])

    def _shape_dxi(self, i):
        return array([-1, 1, 0, 0])

    def _shape_deta(self, i):
        return array([-1, 0, 1, 0])

    def _shape_dpsi(self, i):
        return array([-1, 0, 0, 1])

    # Формирование локальной ыматрицы жесткости
    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((3 * self.size, 3 * self.size))
        self.load = zeros(3 * self.size)
        if not is_static:
            self.M = zeros((3 * self.size, 3 * self.size))
            self.C = zeros((3 * self.size, 3 * self.size))
        # Интегрирование по тетраэдру [0; 0; 0] - [1; 0; 0] - [0; 1; 0] - [0; 0; 1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Матрица Якоби
            jacobi = array([
                [sum(self._shape_dxi(i) * self.x[:, 0]), sum(self._shape_dxi(i) * self.x[:, 1]),
                 sum(self._shape_dxi(i) * self.x[:, 2])],
                [sum(self._shape_deta(i) * self.x[:, 0]), sum(self._shape_deta(i) * self.x[:, 1]),
                 sum(self._shape_deta(i) * self.x[:, 2])],
                [sum(self._shape_dpsi(i) * self.x[:, 0]), sum(self._shape_dpsi(i) * self.x[:, 1]),
                 sum(self._shape_dpsi(i) * self.x[:, 2])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inv_jacobi = inv(jacobi)
            shape_dx = (inv_jacobi[0, 0] * self._shape_dxi(i) + inv_jacobi[0, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[0, 2] * self._shape_dpsi(i))
            shape_dy = (inv_jacobi[1, 0] * self._shape_dxi(i) + inv_jacobi[1, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[1, 2] * self._shape_dpsi(i))
            shape_dz = (inv_jacobi[2, 0] * self._shape_dxi(i) + inv_jacobi[2, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[2, 2] * self._shape_dpsi(i))
            # Изопараметрические матрицы градиентов и функций форм
            b = zeros([6, 3 * self.size])
            n = zeros([3, 3 * self.size])
            for j in range(0, self.size):
                b[0][j * 3 + 0] = b[3][j * 3 + 1] = b[5][j * 3 + 2] = shape_dx[j]
                b[1][j * 3 + 1] = b[3][j * 3 + 0] = b[4][j * 3 + 2] = shape_dy[j]
                b[2][j * 3 + 2] = b[4][j * 3 + 1] = b[5][j * 3 + 0] = shape_dz[j]
                n[0][j * 3 + 0] = n[1][j * 3 + 1] = n[2][j * 3 + 2] = self._shape(i)[j]
            self.K += b.conj().transpose().dot(self._elastic_matrix()).dot(b) * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(v_load) * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += (self._shape(i).conj().transpose().dot(self._shape(i)) * abs(jacobian) * self._w[i] *
                           self.params.density)
                self.C += (self._shape(i).conj().transpose().dot(self._shape(i)) * abs(jacobian) * self._w[i] *
                           self.params.damping)


# Квадратичный (десятиузловой) тетраэдральный КЭ
class TFE3D10(TFE3D4):
    def __init__(self):
        super().__init__()
        self.size = 10

    def _create(self):
        if self._volume() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][2]
                a[i][4] = self.x[i][0] * self.x[i][1]
                a[i][5] = self.x[i][0] * self.x[i][2]
                a[i][6] = self.x[i][1] * self.x[i][2]
                a[i][7] = self.x[i][0] ** 2
                a[i][8] = self.x[i][1] ** 2
                a[i][9] = self.x[i][2] ** 2
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([6, 3 * self.size])
            for j in range(0, self.size):
                dx = (self.a[j][1] + 2 * self.a[j][7] * self.x[i][0] + self.a[j][4] * self.x[i][1] +
                      self.a[j][5] * self.x[i][2])
                dy = (self.a[j][2] + self.a[j][4] * self.x[i][0] + 2 * self.a[j][8] * self.x[i][1] +
                      self.a[j][6] * self.x[i][2])
                dz = (self.a[j][3] + self.a[j][5] * self.x[i][0] + self.a[j][6] * self.x[i][1] +
                      2 * self.a[j][9] * self.x[i][2])
                b[0][j * 3 + 0] = b[3][j * 3 + 1] = b[5][j * 3 + 2] = dx
                b[1][j * 3 + 1] = b[3][j * 3 + 0] = b[4][j * 3 + 2] = dy
                b[2][j * 3 + 2] = b[4][j * 3 + 1] = b[5][j * 3 + 0] = dz

            e = b.dot(u)
            s = self._elastic_matrix().dot(e)
            for j in range(0, 6):
                res[j][i] += e[j]
                res[j + 6][i] += s[j]
        return res

    # Изопараметрические функции формы и их производные
    def _shape(self, i):
        s = array([1 - self._xi[i] - self._eta[i] - self._psi[i], self._xi[i], self._eta[i], self._psi[i]])
        return array([s[0] * (2 * s[0] - 1), s[1] * (2 * s[1] - 1), s[2] * (2 * s[2] - 1), s[3] * (2 * s[3] - 1),
                      4 * s[0] * s[1], 4 * s[1] * s[2], 4 * s[0] * s[2], 4 * s[2] * s[3], 4 * s[1] * s[3],
                      4 * s[0] * s[3]])

    def _shape_dxi(self, i):
        return array([-3 + 4 * self._xi[i] + 4 * self._eta[i] + 4 * self._psi[i], 4 * self._xi[i] - 1, 0, 0,
                      -8 * self._xi[i] + 4 - 4 * self._eta[i] - 4 * self._psi[i],
                      4 * self._eta[i], -4 * self._eta[i], 0, 4 * self._psi[i], -4 * self._psi[i]])

    def _shape_deta(self, i):
        return array([-3 + 4 * self._xi[i] + 4 * self._eta[i] + 4 * self._psi[i], 0, 4 * self._eta[i] - 1, 0,
                      -4 * self._xi[i], 4 * self._xi[i], -8 * self._eta[i] + 4 - 4 * self._xi[i] - 4 * self._psi[i],
                      4 * self._psi[i], 0, -4 * self._psi[i]])

    def _shape_dpsi(self, i):
        return array([-3 + 4 * self._xi[i] + 4 * self._eta[i] + 4 * self._psi[i], 0, 0, 4 * self._psi[i] - 1,
                      -4 * self._xi[i], 0, -4 * self._eta[i], 4 * self._eta[i],
                      4 * self._xi[i], -8 * self._psi[i] + 4 - 4 * self._xi[i] - 4 * self._eta[i]])

    # Формирование локальной ыматрицы жесткости
    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((3 * self.size, 3 * self.size))
        self.load = zeros(3 * self.size)
        if not is_static:
            self.M = zeros((3 * self.size, 3 * self.size))
            self.C = zeros((3 * self.size, 3 * self.size))
        # Интегрирование по тетраэдру [0; 0; 0] - [1; 0; 0] - [0; 1; 0] - [0; 0; 1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Матрица Якоби
            jacobi = array([
                [sum(self._shape_dxi(i) * self.x[:, 0]), sum(self._shape_dxi(i) * self.x[:, 1]),
                 sum(self._shape_dxi(i) * self.x[:, 2])],
                [sum(self._shape_deta(i) * self.x[:, 0]), sum(self._shape_deta(i) * self.x[:, 1]),
                 sum(self._shape_deta(i) * self.x[:, 2])],
                [sum(self._shape_dpsi(i) * self.x[:, 0]), sum(self._shape_dpsi(i) * self.x[:, 1]),
                 sum(self._shape_dpsi(i) * self.x[:, 2])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            if jacobian == 0:
                print('Wrong jacobian: %f' % jacobian)
            inv_jacobi = inv(jacobi)
            shape_dx = (inv_jacobi[0, 0] * self._shape_dxi(i) + inv_jacobi[0, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[0, 2] * self._shape_dpsi(i))
            shape_dy = (inv_jacobi[1, 0] * self._shape_dxi(i) + inv_jacobi[1, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[1, 2] * self._shape_dpsi(i))
            shape_dz = (inv_jacobi[2, 0] * self._shape_dxi(i) + inv_jacobi[2, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[2, 2] * self._shape_dpsi(i))
            # Изопараметрическая матрица градиентов
            b = zeros([6, 3 * self.size])
            # Матрица функций форм
            n = zeros([3, 3 * self.size])
            for j in range(0, self.size):
                b[0][j * 3 + 0] = b[3][j * 3 + 1] = b[5][j * 3 + 2] = shape_dx[j]
                b[1][j * 3 + 1] = b[3][j * 3 + 0] = b[4][j * 3 + 2] = shape_dy[j]
                b[2][j * 3 + 2] = b[4][j * 3 + 1] = b[5][j * 3 + 0] = shape_dz[j]
                n[0][j * 3 + 0] = n[1][j * 3 + 1] = n[2][j * 3 + 2] = self._shape(i)[j]

            self.K += b.conj().transpose().dot(self._elastic_matrix()).dot(b) * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(v_load) * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += (self._shape(i).conj().transpose().dot(self._shape(i)) * abs(jacobian) * self._w[i] *
                           self.params.density)
                self.C += (self._shape(i).conj().transpose().dot(self._shape(i)) * abs(jacobian) * self._w[i] *
                           self.params.damping)


# Восьмиузловой призматический КЭ
class TFE3D8(TFE):
    def __init__(self):
        super().__init__()
        self.size = 8
        # Параметры квадратур Гаусса
        self._xi = [-0.57735027, -0.57735027, -0.57735027, -0.57735027, 0.57735027, 0.57735027, 0.57735027, 0.57735027]
        self._eta = [-0.57735027, -0.57735027, 0.57735027, 0.57735027, -0.57735027, -0.57735027, 0.57735027, 0.57735027]
        self._psi = [-0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027, 0.57735027]
        self._w = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    def _create(self):
        self.a, a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0] * self.size)
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][2]
                a[i][4] = self.x[i][0] * self.x[i][1]
                a[i][5] = self.x[i][0] * self.x[i][2]
                a[i][6] = self.x[i][1] * self.x[i][2]
                a[i][7] = self.x[i][0] * self.x[i][1] * self.x[i][2]
            b[j] = 1.0
            try:
                x = solve(a, b)
            except LinAlgError:
                raise TFEMException('incorrect_fe_err')
            self.a[j] = list(x)

    def _elastic_matrix(self):
        # Матрица упругих свойст
        return array([
            [1.0, self.params.m[0]/(1.0 - self.params.m[0]), self.params.m[0]/(1.0 - self.params.m[0]), 0.0, 0.0, 0.0],
            [self.params.m[0]/(1.0 - self.params.m[0]), 1.0, self.params.m[0]/(1.0 - self.params.m[0]), 0.0, 0.0, 0.0],
            [self.params.m[0]/(1.0 - self.params.m[0]), self.params.m[0]/(1.0 - self.params.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * self.params.m[0])/(1.0 - self.params.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * self.params.m[0])/(1.0 - self.params.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5 * (1.0 - 2.0 * self.params.m[0])/(1.0 - self.params.m[0])],
        ]) * self.params.e[0] * (1.0 - self.params.m[0])/(1.0 + self.params.m[0])/(1.0 - 2.0 * self.params.m[0])

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([6, 3 * self.size])
            for j in range(0, self.size):
                dx = self.a[j][1] + self.a[j][4] * self.x[i][1] + self.a[j][5] * self.x[i][2] + \
                    self.a[j][7] * self.x[i][1] * self.x[i][2]
                dy = self.a[j][2] + self.a[j][4] * self.x[i][0] + self.a[j][6] * self.x[i][2] + \
                    self.a[j][7] * self.x[i][0] * self.x[i][2]
                dz = self.a[j][3] + self.a[j][5] * self.x[i][0] + self.a[j][6] * self.x[i][1] + \
                    self.a[j][7] * self.x[i][0] * self.x[i][1]
                b[0][j * 3 + 0] = b[3][j * 3 + 1] = b[5][j * 3 + 2] = dx
                b[1][j * 3 + 1] = b[3][j * 3 + 0] = b[4][j * 3 + 2] = dy
                b[2][j * 3 + 2] = b[4][j * 3 + 1] = b[5][j * 3 + 0] = dz

            e = b.dot(u)
            s = self._elastic_matrix().dot(e)
            for j in range(0, 6):
                res[j][i] += e[j]
                res[j + 6][i] += s[j]
        return res

    # Изопараметрические функции формы и их производные
    def _shape(self, i):
        return array([
            0.125 * (1.0 - self._xi[i]) * (1.0 - self._eta[i]) * (1.0 - self._psi[i]),
            0.125 * (1.0 + self._xi[i]) * (1.0 - self._eta[i]) * (1.0 - self._psi[i]),
            0.125 * (1.0 + self._xi[i]) * (1.0 + self._eta[i]) * (1.0 - self._psi[i]),
            0.125 * (1.0 - self._xi[i]) * (1.0 + self._eta[i]) * (1.0 - self._psi[i]),
            0.125 * (1.0 - self._xi[i]) * (1.0 - self._eta[i]) * (1.0 + self._psi[i]),
            0.125 * (1.0 + self._xi[i]) * (1.0 - self._eta[i]) * (1.0 + self._psi[i]),
            0.125 * (1.0 + self._xi[i]) * (1.0 + self._eta[i]) * (1.0 + self._psi[i]),
            0.125 * (1.0 - self._xi[i]) * (1.0 + self._eta[i]) * (1.0 + self._psi[i])
        ])

    def _shape_dxi(self, i):
        return array([
            -0.125 * (1.0 - self._eta[i]) * (1.0 - self._psi[i]),
            0.125 * (1.0 - self._eta[i]) * (1.0 - self._psi[i]),
            0.125 * (1.0 + self._eta[i]) * (1.0 - self._psi[i]),
            -0.125 * (1.0 + self._eta[i]) * (1.0 - self._psi[i]),
            -0.125 * (1.0 - self._eta[i]) * (1.0 + self._psi[i]),
            0.125 * (1.0 - self._eta[i]) * (1.0 + self._psi[i]),
            0.125 * (1.0 + self._eta[i]) * (1.0 + self._psi[i]),
            -0.125 * (1.0 + self._eta[i]) * (1.0 + self._psi[i])
        ])

    def _shape_deta(self, i):
        return array([
            -0.125 * (1.0 - self._xi[i]) * (1.0 - self._psi[i]),
            -0.125 * (1.0 + self._xi[i]) * (1.0 - self._psi[i]),
            0.125 * (1.0 + self._xi[i]) * (1.0 - self._psi[i]),
            0.125 * (1.0 - self._xi[i]) * (1.0 - self._psi[i]),
            -0.125 * (1.0 - self._xi[i]) * (1.0 + self._psi[i]),
            -0.125 * (1.0 + self._xi[i]) * (1.0 + self._psi[i]),
            0.125 * (1.0 + self._xi[i]) * (1.0 + self._psi[i]),
            0.125 * (1.0 - self._xi[i]) * (1.0 + self._psi[i])
        ])

    def _shape_dpsi(self, i):
        return array([
            -0.125 * (1.0 - self._xi[i]) * (1.0 - self._eta[i]),
            -0.125 * (1.0 + self._xi[i]) * (1.0 - self._eta[i]),
            -0.125 * (1.0 + self._xi[i]) * (1.0 + self._eta[i]),
            -0.125 * (1.0 - self._xi[i]) * (1.0 + self._eta[i]),
            0.125 * (1.0 - self._xi[i]) * (1.0 - self._eta[i]),
            0.125 * (1.0 + self._xi[i]) * (1.0 - self._eta[i]),
            0.125 * (1.0 + self._xi[i]) * (1.0 + self._eta[i]),
            0.125 * (1.0 - self._xi[i]) * (1.0 + self._eta[i])
        ])

    # Формирование локальной матрицы жесткости
    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((3 * self.size, 3 * self.size))
        self.load = zeros(3 * self.size)
        if not is_static:
            self.M = zeros((3 * self.size, 3 * self.size))
            self.C = zeros((3 * self.size, 3 * self.size))
        # Интегрирование по кубу [-1; 1] x [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Матрица Якоби
            jacobi = array([
                [sum(self._shape_dxi(i) * self.x[:, 0]), sum(self._shape_dxi(i) * self.x[:, 1]),
                 sum(self._shape_dxi(i) * self.x[:, 2])],
                [sum(self._shape_deta(i) * self.x[:, 0]), sum(self._shape_deta(i) * self.x[:, 1]),
                 sum(self._shape_deta(i) * self.x[:, 2])],
                [sum(self._shape_dpsi(i) * self.x[:, 0]), sum(self._shape_dpsi(i) * self.x[:, 1]),
                 sum(self._shape_dpsi(i) * self.x[:, 2])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inv_jacobi = inv(jacobi)
            shape_dx = (inv_jacobi[0, 0] * self._shape_dxi(i) + inv_jacobi[0, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[0, 2] * self._shape_dpsi(i))
            shape_dy = (inv_jacobi[1, 0] * self._shape_dxi(i) + inv_jacobi[1, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[1, 2] * self._shape_dpsi(i))
            shape_dz = (inv_jacobi[2, 0] * self._shape_dxi(i) + inv_jacobi[2, 1] * self._shape_deta(i)) + \
                       (inv_jacobi[2, 2] * self._shape_dpsi(i))
            # Изопараметрическая матрица градиентов
            b = zeros([6, 3 * self.size])
            # Матрица функций форм
            n = zeros([3, 3 * self.size])
            for j in range(0, self.size):
                b[0][j * 3 + 0] = b[3][j * 3 + 1] = b[5][j * 3 + 2] = shape_dx[j]
                b[1][j * 3 + 1] = b[3][j * 3 + 0] = b[4][j * 3 + 2] = shape_dy[j]
                b[2][j * 3 + 2] = b[4][j * 3 + 1] = b[5][j * 3 + 0] = shape_dz[j]
                n[0][j * 3 + 0] = n[1][j * 3 + 1] = n[2][j * 3 + 2] = self._shape(i)[j]
            self.K += b.conj().transpose().dot(self._elastic_matrix()).dot(b) * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(v_load) * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += (self._shape(i).conj().transpose().dot(self._shape(i)) * abs(jacobian) * self._w[i] *
                           self.params.density)
                self.C += (self._shape(i).conj().transpose().dot(self._shape(i)) * abs(jacobian) * self._w[i] *
                           self.params.damping)


# Четырехугольный КЭ пластины
class TFE2D4P(TFE2D4):
    def __init__(self):
        super().__init__()

    def _extra_elastic_matrix(self):
        return array([
            [1.0, 0.0],
            [0.0, 1.0],
        ]) * self.params.e[0]/(2.0 + 2.0 * self.params.m[0])

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b1 = zeros([3, 3 * self.size])
            b2 = zeros([2, 3 * self.size])
            for j in range(0, self.size):
                shape = self.a[j][0] + self.a[j][1] * self.x[i][0] + self.a[j][2] * self.x[i][1] + \
                        self.a[j][3] * self.x[i][0] * self.x[i][1]
                dx = self.a[j][1] + self.a[j][3] * self.x[i][1]
                dy = self.a[j][2] + self.a[j][3] * self.x[i][0]
                b1[0][3 * j + 2] = -dx
                b1[1][3 * j + 1] = dy
                b1[2][3 * j + 1] = dx
                b1[2][3 * j + 2] = -dy
                b2[0][3 * j + 0] = dx
                b2[0][3 * j + 2] = shape
                b2[1][3 * j + 0] = dy
                b2[1][3 * j + 1] = -shape
            e1 = b1.dot(u)
            s1 = self._elastic_matrix().dot(e1)
            e2 = b2.dot(u)
            s2 = self._extra_elastic_matrix().dot(e2)
            res[0][i] += e1[0]  # Exx
            res[1][i] += e1[1]  # Eyy
            res[3][i] += e1[2]  # Exy
            res[4][i] += e2[0]  # Exz
            res[5][i] += e2[1]  # Eyz
            res[6][i] += s1[0]  # Sxx
            res[7][i] += s1[1]  # Syy
            res[9][i] += s1[2]  # Sxy
            res[10][i] += s2[0]  # Sxz
            res[11][i] += s2[1]  # Syz
        return res

    # Формирование локальной матрицы жесткости
    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((3 * self.size, 3 * self.size))
        self.load = zeros(3 * self.size)
        if not is_static:
            self.M = zeros((3 * self.size, 3 * self.size))
            self.C = zeros((3 * self.size, 3 * self.size))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Изопараметрические функции формы и их производные
            shape = self._shape(i)
            shape_dxi = self._shape_dxi(i)
            shape_deta = self._shape_deta(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi * self.x[:, 0]), sum(shape_dxi * self.x[:, 1])],
                [sum(shape_deta * self.x[:, 0]), sum(shape_deta * self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inv_jacobi = inv(jacobi)
            shape_dx = inv_jacobi[0, 0] * shape_dxi + inv_jacobi[0, 1] * shape_deta
            shape_dy = inv_jacobi[1, 0] * shape_dxi + inv_jacobi[1, 1] * shape_deta
            # Изопараметрические матрицы градиентов и функций форм
            b1 = zeros([3, 3 * self.size])
            b2 = zeros([2, 3 * self.size])
            n = zeros([3, 3 * self.size])
            for j in range(0, self.size):
                b1[0][3 * j + 2] = -shape_dx[j]
                b1[1][3 * j + 1] = shape_dy[j]
                b1[2][3 * j + 1] = shape_dx[j]
                b1[2][3 * j + 2] = -shape_dy[j]
                b2[0][3 * j + 0] = shape_dx[j]
                b2[0][3 * j + 2] = shape[j]
                b2[1][3 * j + 0] = shape_dy[j]
                b2[1][3 * j + 1] = -shape[j]
                n[0][3 * j + 0] = n[1][3 * j + 1] = n[2][3 * j + 2] = shape[j]
            self.K += (b1.conj().transpose().dot(self._elastic_matrix()).dot(b1) * self.params.thickness ** 3 / 12.0 +
                       b2.conj().transpose().dot(self._extra_elastic_matrix()).
                       dot(b2) * self.params.thickness * 5.0 / 6.0) * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(v_load) * self.params.thickness * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(s_load) * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += shape.conj().transpose().dot(shape) * abs(jacobian) * self._w[i] * self.params.density
                self.C += shape.conj().transpose().dot(shape) * abs(jacobian) * self._w[i] * self.params.damping


# Треугольный КЭ пластины
class TFE2D3P(TFE2D3):
    def __init__(self):
        super().__init__()

    def _extra_elastic_matrix(self):
        # Вспомогательная матрица упругих свойст
        return array([
            [1.0, 0.0],
            [0.0, 1.0],
        ]) * self.params.e[0]/(2.0 + 2.0 * self.params.m[0])

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b1 = zeros([3, 3 * self.size])
            b2 = zeros([2, 3 * self.size])
            for j in range(0, self.size):
                shape = self.a[j][0] + self.a[j][1] * self.x[i][0] + self.a[j][2] * self.x[i][1]
                dx = self.a[j][1]
                dy = self.a[j][2]
                b1[0][3 * j + 2] = -dx
                b1[1][3 * j + 1] = dy
                b1[2][3 * j + 1] = dx
                b1[2][3 * j + 2] = -dy
                b2[0][3 * j + 0] = dx
                b2[0][3 * j + 2] = shape
                b2[1][3 * j + 0] = dy
                b2[1][3 * j + 1] = -shape
            e1 = b1.dot(u)
            s1 = self._elastic_matrix().dot(e1)
            e2 = b2.dot(u)
            s2 = self._extra_elastic_matrix().dot(e2)
            res[0][i] += e1[0]  # Exx
            res[1][i] += e1[1]  # Eyy
            res[3][i] += e1[2]  # Exy
            res[4][i] += e2[0]  # Exz
            res[5][i] += e2[1]  # Eyz

            res[6][i] += s1[0]  # Sxx
            res[7][i] += s1[1]  # Syy
            res[9][i] += s1[2]  # Sxy
            res[10][i] += s2[0]  # Sxz
            res[11][i] += s2[1]  # Syz
        return res

    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((3 * self.size, 3 * self.size))
        self.load = zeros(3 * self.size)
        if not is_static:
            self.M = zeros((3 * self.size, 3 * self.size))
            self.C = zeros((3 * self.size, 3 * self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[0,1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Изопараметрические функции формы и их производные
            shape = self._shape(i)
            shape_dxi = self._shape_dxi(i)
            shape_deta = self._shape_deta(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi * self.x[:, 0]), sum(shape_dxi * self.x[:, 1])],
                [sum(shape_deta * self.x[:, 0]), sum(shape_deta * self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inv_jacobi = inv(jacobi)
            shape_dx = inv_jacobi[0, 0] * shape_dxi + inv_jacobi[0, 1] * shape_deta
            shape_dy = inv_jacobi[1, 0] * shape_dxi + inv_jacobi[1, 1] * shape_deta
            # Матрицы градиентов и функций форм
            b1 = zeros([3, 3 * self.size])
            b2 = zeros([2, 3 * self.size])
            n = zeros([3, 3 * self.size])
            for j in range(self.size):
                b1[0][3 * j + 2] = -shape_dx[j]
                b1[1][3 * j + 1] = shape_dy[j]
                b1[2][3 * j + 1] = shape_dx[j]
                b1[2][3 * j + 2] = -shape_dy[j]
                b2[0][3 * j + 0] = shape_dx[j]
                b2[0][3 * j + 2] = shape[j]
                b2[1][3 * j + 0] = shape_dy[j]
                b2[1][3 * j + 1] = -shape[j]
                n[0][3 * j + 0] = n[1][3 * j + 1] = n[2][3 * j + 2] = shape[j]
            # Вычисление компонент локальной матрицы жесткости
            self.K += (b1.conj().transpose().dot(self._elastic_matrix()).dot(b1) * self.params.thickness ** 3 / 12.0 +
                       b2.conj().transpose().dot(self._extra_elastic_matrix()).dot(b2) * self.params.thickness *
                       5.0 / 6.0) * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(v_load) * self.params.thickness * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(s_load) * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += shape.conj().transpose().dot(shape) * abs(jacobian) * self._w[i] * self.params.density
                self.C += shape.conj().transpose().dot(shape) * abs(jacobian) * self._w[i] * self.params.damping


# Треугольный КЭ оболочки
class TFE2D3S(TFE2D3P):
    def __init__(self):
        super().__init__()
        self.size = 3
        self.T = zeros((3, 3))  # Матрица преобразования глобальных координат в локальные
        self.global_x = zeros((3, 3))

    def _create(self):
        self.T = create_transform_matrix(self.x)
        self.global_x = self.x
        self.x = array([self.T.dot(self.x[0, :]), self.T.dot(self.x[1, :]), self.T.dot(self.x[2, :])]) - self.x[0, :]
        TFE2D3._create(self)

    def _prepare_transform_matrix_(self):
        m = zeros((6 * self.size, 6 * self.size))
        for i in range(0, 3):
            for j in range(0, 3):
                m[i][j] = self.T[i][j]
                m[i + 3][j + 3] = self.T[i][j]
                m[i + 6][j + 6] = self.T[i][j]
                m[i + 9][j + 9] = self.T[i][j]
                m[i + 12][j + 12] = self.T[i][j]
                m[i + 15][j + 15] = self.T[i][j]
        return m

    def calc(self, u):
        # Подготовка матрицы преобразования
        m = self._prepare_transform_matrix_()
        # Модифицируем перемещения
        lu = m.dot(u)
        # Вычисляем деформации и напряжения "мембранной" составляющей
        membrane_u = []
        for i in range(0, self.size):
            membrane_u.append(lu[6 * i])      # u
            membrane_u.append(lu[6 * i + 1])  # v
        membrane_res = TFE2D3.calc(self, membrane_u)
        # Вычисляем деформации и напряжения пластины
        plate_u = []
        for i in range(0, self.size):
            plate_u.append(lu[6 * i + 2])     # w
            plate_u.append(lu[6 * i + 3])     # Tx
            plate_u.append(lu[6 * i + 4])     # Ty
        plate_res = TFE2D3P.calc(self, plate_u)
        res = zeros((12, self.size))
        for i in range(0, self.size):
            ld = array([
                [plate_res[0][i] + membrane_res[0][i], plate_res[3][i] + membrane_res[2][i], plate_res[4][i]],
                [plate_res[3][i] + membrane_res[2][i], plate_res[1][i] + membrane_res[1][i], plate_res[5][i]],
                [plate_res[4][i], plate_res[5][i], 0]
            ])
            ls = array([
                [plate_res[6][i] + membrane_res[3][i], plate_res[9][i] + membrane_res[5][i], plate_res[10][i]],
                [plate_res[9][i] + membrane_res[5][i], plate_res[7][i] + membrane_res[4][i], plate_res[11][i]],
                [plate_res[10][i], plate_res[11][i], 0]
            ])
            gd = self.T.conj().transpose().dot(ld).dot(self.T)
            gs = self.T.conj().transpose().dot(ls).dot(self.T)

            res[0][i] = gd[0][0]    # Exx
            res[1][i] = gd[1][1]    # Eyy
            res[2][i] = gd[2][2]    # Ezz
            res[3][i] = gd[0][1]    # Exy
            res[4][i] = gd[0][2]    # Exz
            res[5][i] = gd[1][2]    # Eyz
            res[6][i] = gs[0][0]    # Sxx
            res[7][i] = gs[1][1]    # Syy
            res[8][i] = gs[2][2]    # Szz
            res[9][i] = gs[0][1]    # Sxy
            res[10][i] = gs[0][2]   # Sxz
            res[11][i] = gs[1][2]   # Syz
        return res

    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((6 * self.size, 6 * self.size))
        self.load = zeros(6 * self.size)
        if not is_static:
            self.M = zeros((6 * self.size, 6 * self.size))
            self.C = zeros((6 * self.size, 6 * self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[0,1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Изопараметрические функции формы и их производные
            shape = self._shape(i)
            shape_dxi = self._shape_dxi(i)
            shape_deta = self._shape_deta(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi * self.x[:, 0]), sum(shape_dxi * self.x[:, 1])],
                [sum(shape_deta * self.x[:, 0]), sum(shape_deta * self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inv_jacobi = inv(jacobi)
            shape_dx = inv_jacobi[0, 0] * shape_dxi + inv_jacobi[0, 1] * shape_deta
            shape_dy = inv_jacobi[1, 0] * shape_dxi + inv_jacobi[1, 1] * shape_deta
            # Матрицы градиентов и функций форм
            bm = zeros([3, 6 * self.size])
            bb = zeros([3, 6 * self.size])
            bc = zeros([2, 6 * self.size])
            n = zeros([6, 6 * self.size])
            for j in range(self.size):
                bm[0][6 * j + 0] = shape_dx[j]
                bm[1][6 * j + 1] = shape_dy[j]
                bm[2][6 * j + 0] = shape_dy[j]
                bm[2][6 * j + 1] = shape_dx[j]
                bb[0][6 * j + 3] = shape_dx[j]
                bb[1][6 * j + 4] = shape_dy[j]
                bb[2][6 * j + 3] = shape_dy[j]
                bb[2][6 * j + 4] = shape_dx[j]
                bc[0][6 * j + 2] = shape_dx[j]
                bc[0][6 * j + 3] = shape[j]
                bc[1][6 * j + 0] = shape_dy[j]
                bc[1][6 * j + 4] = shape[j]

                n[0][6 * j + 0] = n[1][6 * j + 1] = n[2][6 * j + 2] = n[3][6 * j + 3] = n[4][6 * j + 4] = \
                    n[5][6 * j + 5] = shape[j]

            # Вычисление компонент локальной матрицы жесткости
            self.K += (bm.conj().transpose().dot(self._elastic_matrix()).dot(bm) * self.params.thickness +
                       bb.conj().transpose().dot(self._elastic_matrix()).dot(bb) * self.params.thickness ** 3 / 12.0 +
                       bc.conj().transpose().dot(self._extra_elastic_matrix()).dot(bc) * self.params.thickness *
                       5 / 6) * abs(jacobian) * self._w[i]
            # Вычисление столбца нагрузки
            self.load += n.conj().transpose().dot(v_load) * self.params.thickness * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(s_load) * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += shape.conj().transpose().dot(shape) * abs(jacobian) * self._w[i] * self.params.density
                self.C += shape.conj().transpose().dot(shape) * abs(jacobian) * self._w[i] * self.params.damping

        # Устранение сингулярности
        for i in range(len(self.K)):
            if self.K[i][i] == 0:
                self.K[i][i] = 1.0E+10
                if not is_static:
                    self.M[i][i] = self.C[i][i] = 1.0E+10

        # Подготовка матрицы преобразования
        m = self._prepare_transform_matrix_()
        self.K = m.conj().transpose().dot(self.K).dot(m)
        if not is_static:
            self.M = m.conj().transpose().dot(self.M).dot(m)
            self.C = m.conj().transpose().dot(self.C).dot(m)


# Четырехугольный КЭ оболочки
class TFE2D4S(TFE2D4P):
    def __init__(self):
        super().__init__()
        self.T = zeros((3, 3))  # Матрица преобразования глобальных координат в локальные
        self.global_x = zeros((4, 4))

    def calc(self, u):
        # Подготовка матрицы преобразования
        m = self._prepare_transform_matrix_()
        # Модифицируем перемещения
        lu = m.dot(u)
        # Вычисляем деформации и напряжения "мембранной" составляющей
        membrane_u = []
        for i in range(0, self.size):
            membrane_u.append(lu[6 * i])      # u
            membrane_u.append(lu[6 * i + 1])  # w
        membrane_res = TFE2D4.calc(self, membrane_u)
        # Вычисляем деформации и напряжения пластины
        plate_u = []
        for i in range(0, self.size):
            plate_u.append(lu[6 * i + 2])     # w
            plate_u.append(lu[6 * i + 3])     # Tx
            plate_u.append(lu[6 * i + 4])     # Ty
        plate_res = TFE2D4P.calc(self, plate_u)
        res = zeros((12, self.size))
        for i in range(0, self.size):
            ld = array([
                [plate_res[0][i] + membrane_res[0][i], plate_res[3][i] + membrane_res[2][i], plate_res[4][i]],
                [plate_res[3][i] + membrane_res[2][i], plate_res[1][i] + membrane_res[1][i], plate_res[5][i]],
                [plate_res[4][i], plate_res[5][i], 0]
            ])
            ls = array([
                [plate_res[6][i] + membrane_res[3][i], plate_res[9][i] + membrane_res[5][i], plate_res[10][i]],
                [plate_res[9][i] + membrane_res[5][i], plate_res[7][i] + membrane_res[4][i], plate_res[11][i]],
                [plate_res[10][i], plate_res[11][i], 0]
            ])
            gd = self.T.conj().transpose().dot(ld).dot(self.T)
            gs = self.T.conj().transpose().dot(ls).dot(self.T)

            res[0][i] = gd[0][0]    # Exx
            res[1][i] = gd[1][1]    # Eyy
            res[2][i] = gd[2][2]    # Ezz
            res[3][i] = gd[0][1]    # Exy
            res[4][i] = gd[0][2]    # Exz
            res[5][i] = gd[1][2]    # Eyz
            res[6][i] = gs[0][0]    # Sxx
            res[7][i] = gs[1][1]    # Syy
            res[8][i] = gs[2][2]    # Szz
            res[9][i] = gs[0][1]    # Sxy
            res[10][i] = gs[0][2]   # Sxz
            res[11][i] = gs[1][2]   # Syz
        return res

    def _prepare_transform_matrix_(self):
        m = zeros((6 * self.size, 6 * self.size))
        for i in range(0, 3):
            for j in range(0, 3):
                m[i][j] = self.T[i][j]
                m[i + 3][j + 3] = self.T[i][j]
                m[i + 6][j + 6] = self.T[i][j]
                m[i + 9][j + 9] = self.T[i][j]
                m[i + 12][j + 12] = self.T[i][j]
                m[i + 15][j + 15] = self.T[i][j]
                m[i + 18][j + 18] = self.T[i][j]
                m[i + 21][j + 21] = self.T[i][j]
        return m

    # Формирование локальной матрицы жесткости
    def generate(self, v_load, s_load, is_static=True):
        self.K = zeros((6 * self.size, 6 * self.size))
        self.load = zeros(6 * self.size)
        if not is_static:
            self.M = zeros((6 * self.size, 6 * self.size))
            self.C = zeros((6 * self.size, 6 * self.size))

        # Интегрирование по квадрату [-1,-1]-[1,1] (по формуле Гаусса)
        for i in range(len(self._w)):
            # Изопараметрические функции формы и их производные
            shape = self._shape(i)
            shape_dxi = self._shape_dxi(i)
            shape_deta = self._shape_deta(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi * self.x[:, 0]), sum(shape_dxi * self.x[:, 1])],
                [sum(shape_deta * self.x[:, 0]), sum(shape_deta * self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inv_jacobi = inv(jacobi)
            shape_dx = inv_jacobi[0, 0] * shape_dxi + inv_jacobi[0, 1] * shape_deta
            shape_dy = inv_jacobi[1, 0] * shape_dxi + inv_jacobi[1, 1] * shape_deta
            # Матрицы градиентов и функций форм
            bm = zeros([3, 6 * self.size])
            bb = zeros([3, 6 * self.size])
            bc = zeros([2, 6 * self.size])
            n = zeros([6, 6 * self.size])
            for j in range(self.size):
                bm[0][6 * j + 0] = shape_dx[j]
                bm[1][6 * j + 1] = shape_dy[j]
                bm[2][6 * j + 0] = shape_dy[j]
                bm[2][6 * j + 1] = shape_dx[j]

                bb[0][6 * j + 3] = shape_dy[j]
                bb[1][6 * j + 4] = shape_dy[j]
                bb[2][6 * j + 3] = shape_dy[j]
                bb[2][6 * j + 4] = shape_dx[j]

                bc[0][6 * j + 2] = shape_dx[j]
                bc[0][6 * j + 3] = shape[j]
                bc[1][6 * j + 2] = shape_dy[j]
                bc[1][6 * j + 4] = shape[j]
                n[0][6 * j + 0] = n[1][6 * j + 1] = n[2][6 * j + 2] = n[3][6 * j + 3] = n[4][6 * j + 4] = \
                    n[5][6 * j + 5] = shape[j]

            # Вычисление компонент локальной матрицы жесткости
            self.K += (bm.conj().transpose().dot(self._elastic_matrix()).dot(bm) * self.params.thickness +
                       bb.conj().transpose().dot(self._elastic_matrix()).dot(bb) * (self.params.thickness ** 3) / 12.0 +
                       bc.conj().transpose().dot(self._extra_elastic_matrix()).dot(bc) * self.params.thickness *
                       5.0 / 6.0) * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(v_load) * self.params.thickness * abs(jacobian) * self._w[i]
            self.load += n.conj().transpose().dot(s_load) * abs(jacobian) * self._w[i]
            if not is_static:
                self.M += shape.conj().transpose().dot(shape) * abs(jacobian) * self._w[i] * self.params.density
                self.C += shape.conj().transpose().dot(shape) * abs(jacobian) * self._w[i] * self.params.damping

        # Устранение сингулярности матрицы
        for i in range(len(self.K)):
            if self.K[i][i] == 0:
                self.K[i][i] = 1.0E+10

        # Создание матрицы преобразования
        m = self._prepare_transform_matrix_()
        self.K = m.conj().transpose().dot(self.K).dot(m)


    def _create(self):
        self.T = create_transform_matrix(self.x)
        self.global_x = self.x

        x0 = self.T.dot(self.x[0, :].conj().transpose())
        x1 = self.T.dot(self.x[1, :].conj().transpose())
        x2 = self.T.dot(self.x[2, :].conj().transpose())
        x3 = self.T.dot(self.x[3, :].conj().transpose())

        self.x = array([x0 - x0, x1 - x0, x2 - x0, x3 - x0])
#        self.x = array([x0, x1, x2, x3])
        TFE2D4._create(self)
