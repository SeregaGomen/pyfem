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


# Абстрактный базовый класс, описывающий конечный элемент (КЭ)
class TFE:
    def __init__(self):
        self.size = 0           # Размерность КЭ
        self.h = 1              # Толщина (для оболочек и пластин)
        self.density = 0        # Плотность
        self.damping = [0, 0]   # Параметры демпфирования
        self.e = []             # Модуль Юнга
        self.m = []             # Коэффициент Пуассона
        self.x = []             # Координаты вершин КЭ
        self.K = []             # Локальная матрица жесткости
        self.M = []             # ... масс
        self.C = []             # ... демпфирования
        self.a = []             # Коэффициенты функций форм

    # Задание толщины
    def set_h(self, h):
        self.h = h

    # Задание координат
    def set_coord(self, x):
        self.x = array(x)
        self._create_()

    # Задание параметров упругости
    def set_elasticity(self, p1, p2):
        self.e = p1
        self.m = p2

    # Задание плотности
    def set_density(self, d):
        self.density = d

    # Задание параметра демпфирования
    def set_damping(self, d):
        self.damping = d

    # Формирование матриц жесткости, масс и демпфирования
    def generate(self, is_static=True):
        self.K = self._generate_stiffness_matrix_()
        if not is_static:
            self.M, self.C = self._generate_mass_damping_matrix_()

    # Построение вектора для заданной стороны элемента
    def _vector_(self, i, j):
        v = array(self.x[j]) - array(self.x[i])
        return v/norm(v)

    # Матрица преобразования в локальную систему координат
    def _create_transform_matrix_(self):
        v_x = self._vector_(1, 0)
        v_z = self._cross_product_(self._vector_(1, 0), self._vector_(2, 0))
        v_y = self._cross_product_(v_z, v_x)
        return array([v_x, v_y, v_z])

    # Векторное произведение a и b
    @staticmethod
    def _cross_product_(a, b):
        v = cross(a, b)
        return v/norm(v)

    # Вычисление функций форм КЭ
    @abstractmethod
    def _create_(self):
        raise NotImplementedError('Method TFE._create_() is pure virtual')

    # Матрица упругих свойств
    @abstractmethod
    def _elastic_matrix_(self):
        raise NotImplementedError('Method TFE._elastic_matrix_() is pure virtual')

    # Вычисления стандартных результатов КЭ
    @abstractmethod
    def calc(self, u):
        raise NotImplementedError('Method TFE.calc() is pure virtual')

    # Вычисление матрицы жесткости
    @abstractmethod
    def _generate_stiffness_matrix_(self):
        raise NotImplementedError('Method TFE._generate_stiffness_matrix_() is pure virtual')

    # Вычисление матриц масс и демпфирования
    @abstractmethod
    def _generate_mass_damping_matrix_(self):
        raise NotImplementedError('Method TFE._generate_mass_damping_matrix_() is pure virtual')


# Линейный (двухузловой) одномерный КЭ
class TFE1D2(TFE):
    def __init__(self):
        super().__init__()
        self.size = 2

    def _length_(self):
        return math.fabs(self.x[1][0] - self.x[0][0])

    def _create_(self):
        if self._length_() == 0.0:
            raise TFEMException('incorrect_fe_err')
        self.a = zeros((self.size, self.size))
        self.a[0][0] = self.x[1][0]/(self.x[1][0] - self.x[0][0])
        self.a[0][1] = -1.0/(self.x[1][0] - self.x[0][0])
        self.a[1][0] = self.x[0][0]/(self.x[0][0] - self.x[1][0])
        self.a[1][1] = -1.0/(self.x[0][0] - self.x[1][0])

    def calc(self, u):
        res = zeros((self.size, self.size))
        res[0][0] = res[0][1] = u[0]*self.a[0][1] + u[1]*self.a[1][1]
        res[1][0] = res[1][1] = self.e[0]*(u[0]*self.a[0][1] + u[1]*self.a[1][1])
        return res

    def _elastic_matrix_(self):
        return array([
            [1.0, -1.0],
            [-1.0, .0]
        ])*self.e[0]

    def _generate_stiffness_matrix_(self):
        return self._elastic_matrix_()/self._length_()

    # Формирование матрицы демпфирования по Релею
    def _generate_mass_damping_matrix_(self):
        a = array([
            [2.0, 1.0],
            [1.0, 2.0]
        ])
        m = a * self._length_() / 6.0 * self.density
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Линейный (трехузловой) треугольный КЭ
class TFE2D3(TFE):
    def __init__(self):
        super().__init__()
        self.size = 3

    def _elastic_matrix_(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
        ])*self.e[0]/(1.0 - self.m[0]**2)
        return d

    def _gradient_matrix_(self):
        # Матрица градиентов
        b = zeros([3, 2*self.size])
        for i in range(0, self.size):
            # Производные функций формы
            dx = self.a[i][1]
            dy = self.a[i][2]
            b[0][2*i + 0] = dx
            b[1][2*i + 1] = dy
            b[2][2*i + 0] = dy
            b[2][2*i + 1] = dx
        return b

    def _square_(self):
        a = math.sqrt((self.x[0][0] - self.x[1][0])**2 + (self.x[0][1] - self.x[1][1])**2)
        b = math.sqrt((self.x[0][0] - self.x[2][0])**2 + (self.x[0][1] - self.x[2][1])**2)
        c = math.sqrt((self.x[2][0] - self.x[1][0])**2 + (self.x[2][1] - self.x[1][1])**2)
        p = 0.5*(a + b + c)
        return math.sqrt(p*(p - a)*(p - b)*(p - c))

    def _create_(self):
        det0 = self.x[2][1]*self.x[1][0] - self.x[2][1]*self.x[0][0] - self.x[0][1]*self.x[1][0] - \
               self.x[1][1]*self.x[2][0] + self.x[1][1]*self.x[0][0] + self.x[0][1]*self.x[2][0]
        if math.fabs(det0) < eps:
            raise TFEMException('incorrect_fe_err')
        index = [[2, 1], [0, 2], [1, 0]]
        self.a = zeros((self.size, self.size))
        for i in range(0, self.size):
            det1 = self.x[index[i][0]][1]*self.x[index[i][1]][0] - self.x[index[i][1]][1]*self.x[index[i][0]][0]
            det2 = self.x[index[i][1]][1] - self.x[index[i][0]][1]
            det3 = self.x[index[i][0]][0] - self.x[index[i][1]][0]
            self.a[i][0] = det1/det0
            self.a[i][1] = det2/det0
            self.a[i][2] = det3/det0

    def calc(self, u):
        e = self._gradient_matrix_().dot(u)
        s = self._elastic_matrix_().dot(e)
        res = zeros((6, self.size))
        for i in range(0, 3):
            for j in range(0, self.size):
                res[i][j] = e[i]
                res[i + 3][j] = s[i]
        return res

    # Формирование локальной матриц жесткости
    def _generate_stiffness_matrix_(self):
        return self._gradient_matrix_().conj().transpose().dot(self._elastic_matrix_()).\
                   dot(self._gradient_matrix_())*self._square_()*self.h

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        a = array([
            [0.5, 0.0, 0.25, 0.0, 0.25, 0.0],
            [0.0, 0.5, 0.0, 0.25, 0.0, 0.25],
            [0.25, 0.0, 0.5, 0.0, 0.25, 0.0],
            [0.0, 0.25, 0.0, 0.5, 0.0, 0.25],
            [0.25, 0.0, 0.25, 0.0, 0.5, 0.0],
            [0.0, 0.25, 0.0, 0.25, 0.0, 0.5]
        ])
        m = a*self.density*self._square_()
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Линейный (четырехузловой) тетраэдральный КЭ
class TFE3D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4

    def _elastic_matrix_(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), 1.0, self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0])],
        ])*self.e[0]*(1.0 - self.m[0])/(1.0 + self.m[0])/(1.0 - 2.0*self.m[0])
        return d

    def _gradient_matrix_(self):
        # Матрица градиентов
        b = zeros([6, 3*self.size])
        for i in range(0, self.size):
            # Производные функций формы
            dx = self.a[i][1]
            dy = self.a[i][2]
            dz = self.a[i][3]

            b[0][3*i + 0] = dx
            b[1][3*i + 1] = dy
            b[2][3*i + 2] = dz
            b[3][3*i + 0] = dy
            b[3][3*i + 1] = dx
            b[4][3*i + 1] = dz
            b[4][3*i + 2] = dy
            b[5][3*i + 0] = dz
            b[5][3*i + 2] = dx
        return b

    def _volume_(self):
        a = (self.x[1][0] - self.x[0][0])*(self.x[2][1] - self.x[0][1])*(self.x[3][2] - self.x[0][2]) + \
            (self.x[3][0] - self.x[0][0])*(self.x[1][1] - self.x[0][1])*(self.x[2][2] - self.x[0][2]) + \
            (self.x[2][0] - self.x[0][0])*(self.x[3][1] - self.x[0][1])*(self.x[1][2] - self.x[0][2])
        b = (self.x[3][0] - self.x[0][0])*(self.x[2][1] - self.x[0][1])*(self.x[1][2] - self.x[0][2]) + \
            (self.x[2][0] - self.x[0][0])*(self.x[1][1] - self.x[0][1])*(self.x[3][2] - self.x[0][2]) + \
            (self.x[1][0] - self.x[0][0])*(self.x[3][1] - self.x[0][1])*(self.x[2][2] - self.x[0][2])
        return math.fabs(a - b)/6.0

    def _create_(self):
        if self._volume_() == 0.0:
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
        e = self._gradient_matrix_().dot(u)
        s = self._elastic_matrix_().dot(e)
        res = zeros((12, self.size))
        for i in range(0, 6):
            for j in range(0, self.size):
                res[i][j] = e[i]
                res[i + 6][j] = s[i]
        return res

    # Формирование локальной матрицы жесткости
    def _generate_stiffness_matrix_(self):
        return self._gradient_matrix_().conj().transpose().dot(self._elastic_matrix_()).\
                   dot(self._gradient_matrix_())*self._volume_()

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        a = array([
            [1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5],
            [0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5],
            [0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5],
            [0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0]
        ])*0.1
        m = a*self.density*self._volume_()
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Билинейный четырехузловой двумерный КЭ
class TFE2D4(TFE):
    def __init__(self):
        super().__init__()
        self.size = 4
        # Параметры квадратур Гаусса
        self._xi_ = [-0.57735027, -0.57735027, 0.57735027, 0.57735027]
        self._eta_ = [-0.57735027, 0.57735027, -0.57735027, 0.57735027]
        self._w_ = [1.0, 1.0, 1.0, 1.0]

    def _square_(self):
        return math.sqrt((self.x[0][0] - self.x[1][0])**2 + (self.x[0][1] - self.x[1][1])**2)

    def _create_(self):
        if self._square_() == 0.0:
            raise TFEMException('incorrect_fe_err')
        a, self.a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0, 0.0, 0.0, 0.0])
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][0]*self.x[i][1]
            b[j] = 1.0
            x = solve(a, b)
            self.a[j] = list(x)

    def _elastic_matrix_(self):
        # Матрица упругих свойст
        d = array([
            [1.0, self.m[0], 0.0],
            [self.m[0], 1.0, 0.0],
            [0.0, 0.0, 0.5*(1.0 - self.m[0])]
        ])*self.e[0]/(1.0 - self.m[0]**2)
        return d

    def calc(self, u):
        res = zeros((6, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([3, 2*self.size])
            for j in range(0, self.size):
                dx = self.a[j][1] + self.a[j][3] * self.x[i][1]
                dy = self.a[j][2] + self.a[j][3] * self.x[i][0]
                b[0][j*2 + 0] = dx
                b[1][j*2 + 1] = dy
                b[2][j*2 + 0] = dy
                b[2][j*2 + 1] = dx
            e = b.dot(u)
            s = self._elastic_matrix_().dot(e)
            for j in range(0, 3):
                res[j][i] += e[j]
                res[j + 3][i] += s[j]
        return res

    # Изопараметрические функции формы и их производные
    def _shape_(self, i):
        return array([
            0.25*(1.0 - self._xi_[i])*(1.0 - self._eta_[i]),
            0.25*(1.0 + self._xi_[i])*(1.0 - self._eta_[i]),
            0.25*(1.0 + self._xi_[i])*(1.0 + self._eta_[i]),
            0.25*(1.0 - self._xi_[i])*(1.0 + self._eta_[i])
        ])

    def _shape_dxi_(self, i):
        return array([
            -0.25*(1.0 - self._eta_[i]),
            0.25*(1.0 - self._eta_[i]),
            0.25*(1.0 + self._eta_[i]),
            -0.25*(1.0 + self._eta_[i])
        ])

    def _shape_deta_(self, i):
        return array([
            -0.25*(1.0 - self._xi_[i]),
            -0.25*(1.0 + self._xi_[i]),
            0.25*(1.0 + self._xi_[i]),
            0.25*(1.0 - self._xi_[i])
        ])

    # Формирование локальной матрицы жесткости
    def _generate_stiffness_matrix_(self):
        local_k = zeros((8, 8))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w_)):
            # Изопараметрические функции формы и их производные
            shape_dxi = self._shape_dxi_(i)
            shape_deta = self._shape_deta_(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta
            shape_dy = inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta
            # Изопараметрическвая матрица градиентов
            b = zeros([3, 2*self.size])
            for j in range(0, self.size):
                b[0][j*2 + 0] = shape_dx[j]
                b[1][j*2 + 1] = shape_dy[j]
                b[2][j*2 + 0] = shape_dy[j]
                b[2][j*2 + 1] = shape_dx[j]
            local_k += b.conj().transpose().dot(self._elastic_matrix_()).dot(b)*abs(jacobian)*self._w_[i]*self.h
        return local_k

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        local_m = zeros((2*self.size, 2*self.size))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w_)):
            # Изопараметрические функции формы и их производные
            shape = self._shape_(i)
            shape_dxi = self._shape_dxi_(i)
            shape_deta = self._shape_deta_(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                [shape[0], 0.0, shape[1], 0.0, shape[2], 0.0, shape[3], 0.0],
                [0.0, shape[0], 0.0, shape[1], 0.0, shape[2], 0.0, shape[3]]
            ])
            local_m += c.conj().transpose().dot(c)*abs(jacobian)*self._w_[i]
        m = self.density*local_m
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Восьмиузловой призматический КЭ
class TFE3D8(TFE):
    def __init__(self):
        super().__init__()
        self.size = 8
        # Параметры квадратур Гаусса
        self._xi_ = [-0.57735027, -0.57735027, -0.57735027, -0.57735027, 0.57735027, 0.57735027, 0.57735027, 0.57735027]
        self._eta_ = [-0.57735027, -0.57735027, 0.57735027, 0.57735027, -0.57735027, -0.57735027, 0.57735027, 0.57735027]
        self._psi_ = [-0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027, 0.57735027]
        self._w_ = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    def _create_(self):
        self.a, a = zeros((self.size, self.size)), zeros((self.size, self.size))
        for j in range(0, self.size):
            b = array([0.0]*self.size)
            for i in range(0, self.size):
                a[i][0] = 1.0
                a[i][1] = self.x[i][0]
                a[i][2] = self.x[i][1]
                a[i][3] = self.x[i][2]
                a[i][4] = self.x[i][0]*self.x[i][1]
                a[i][5] = self.x[i][0]*self.x[i][2]
                a[i][6] = self.x[i][1]*self.x[i][2]
                a[i][7] = self.x[i][0]*self.x[i][1]*self.x[i][2]
            b[j] = 1.0
            try:
                x = solve(a, b)
            except LinAlgError:
                raise TFEMException('incorrect_fe_err')
            self.a[j] = list(x)

    def _elastic_matrix_(self):
        # Матрица упругих свойст
        return array([
            [1.0, self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), 1.0, self.m[0]/(1.0 - self.m[0]), 0.0, 0.0, 0.0],
            [self.m[0]/(1.0 - self.m[0]), self.m[0]/(1.0 - self.m[0]), 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0]), 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1.0 - 2.0*self.m[0])/(1.0 - self.m[0])],
        ])*self.e[0]*(1.0 - self.m[0])/(1.0 + self.m[0])/(1.0 - 2.0*self.m[0])

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b = zeros([6, 3*self.size])
            for j in range(0, self.size):
                dx = self.a[j][1] + self.a[j][4] * self.x[i][1] + self.a[j][5] * self.x[i][2] + \
                    self.a[j][7] * self.x[i][1] * self.x[i][2]
                dy = self.a[j][2] + self.a[j][4] * self.x[i][0] + self.a[j][6] * self.x[i][2] + \
                    self.a[j][7] * self.x[i][0] * self.x[i][2]
                dz = self.a[j][3] + self.a[j][5] * self.x[i][0] + self.a[j][6] * self.x[i][1] + \
                    self.a[j][7] * self.x[i][0] * self.x[i][1]
                b[0][j*3 + 0] = dx
                b[1][j*3 + 1] = dy
                b[2][j*3 + 2] = dz
                b[3][j*3 + 0] = dy
                b[3][j*3 + 1] = dx
                b[4][j*3 + 1] = dz
                b[4][j*3 + 2] = dy
                b[5][j*3 + 0] = dz
                b[5][j*3 + 2] = dx
            e = b.dot(u)
            s = self._elastic_matrix_().dot(e)
            for j in range(0, 6):
                res[j][i] += e[j]
                res[j + 6][i] += s[j]
        return res

    # Изопараметрические функции формы и их производные
    def _shape_(self, i):
        return array([
            0.125*(1.0 - self._xi_[i])*(1.0 - self._eta_[i])*(1.0 - self._psi_[i]),
            0.125*(1.0 + self._xi_[i])*(1.0 - self._eta_[i])*(1.0 - self._psi_[i]),
            0.125*(1.0 + self._xi_[i])*(1.0 + self._eta_[i])*(1.0 - self._psi_[i]),
            0.125*(1.0 - self._xi_[i])*(1.0 + self._eta_[i])*(1.0 - self._psi_[i]),
            0.125*(1.0 - self._xi_[i])*(1.0 - self._eta_[i])*(1.0 + self._psi_[i]),
            0.125*(1.0 + self._xi_[i])*(1.0 - self._eta_[i])*(1.0 + self._psi_[i]),
            0.125*(1.0 + self._xi_[i])*(1.0 + self._eta_[i])*(1.0 + self._psi_[i]),
            0.125*(1.0 - self._xi_[i])*(1.0 + self._eta_[i])*(1.0 + self._psi_[i])
        ])

    def _shape_dxi_(self, i):
        return array([
            -0.125*(1.0 - self._eta_[i])*(1.0 - self._psi_[i]),
            0.125*(1.0 - self._eta_[i])*(1.0 - self._psi_[i]),
            0.125*(1.0 + self._eta_[i])*(1.0 - self._psi_[i]),
            -0.125*(1.0 + self._eta_[i])*(1.0 - self._psi_[i]),
            -0.125*(1.0 - self._eta_[i])*(1.0 + self._psi_[i]),
            0.125*(1.0 - self._eta_[i])*(1.0 + self._psi_[i]),
            0.125*(1.0 + self._eta_[i])*(1.0 + self._psi_[i]),
            -0.125*(1.0 + self._eta_[i])*(1.0 + self._psi_[i])
        ])

    def _shape_deta_(self, i):
        return array([
            -0.125*(1.0 - self._xi_[i])*(1.0 - self._psi_[i]),
            -0.125*(1.0 + self._xi_[i])*(1.0 - self._psi_[i]),
            0.125*(1.0 + self._xi_[i])*(1.0 - self._psi_[i]),
            0.125*(1.0 - self._xi_[i])*(1.0 - self._psi_[i]),
            -0.125*(1.0 - self._xi_[i])*(1.0 + self._psi_[i]),
            -0.125*(1.0 + self._xi_[i])*(1.0 + self._psi_[i]),
            0.125*(1.0 + self._xi_[i])*(1.0 + self._psi_[i]),
            0.125*(1.0 - self._xi_[i])*(1.0 + self._psi_[i])
        ])

    def _shape_dpsi_(self, i):
        return array([
            -0.125*(1.0 - self._xi_[i])*(1.0 - self._eta_[i]),
            -0.125*(1.0 + self._xi_[i])*(1.0 - self._eta_[i]),
            -0.125*(1.0 + self._xi_[i])*(1.0 + self._eta_[i]),
            -0.125*(1.0 - self._xi_[i])*(1.0 + self._eta_[i]),
            0.125*(1.0 - self._xi_[i])*(1.0 - self._eta_[i]),
            0.125*(1.0 + self._xi_[i])*(1.0 - self._eta_[i]),
            0.125*(1.0 + self._xi_[i])*(1.0 + self._eta_[i]),
            0.125*(1.0 - self._xi_[i])*(1.0 + self._eta_[i])
        ])

    # Формирование локальной ыматрицы жесткости
    def _generate_stiffness_matrix_(self):
        local_k = zeros((3*self.size, 3*self.size))
        # Интегрирование по кубу [-1; 1] x [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w_)):
            # Изопараметрические функции формы и их производные
            shape_dxi = self._shape_dxi_(i)
            shape_deta = self._shape_deta_(i)
            shape_dpsi = self._shape_dpsi_(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1]), sum(shape_dxi*self.x[:, 2])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1]), sum(shape_deta*self.x[:, 2])],
                [sum(shape_dpsi*self.x[:, 0]), sum(shape_dpsi*self.x[:, 1]), sum(shape_dpsi*self.x[:, 2])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = (inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta) + \
                       (inverted_jacobi[0, 2]*shape_dpsi)
            shape_dy = (inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta) + \
                       (inverted_jacobi[1, 2]*shape_dpsi)
            shape_dz = (inverted_jacobi[2, 0]*shape_dxi + inverted_jacobi[2, 1]*shape_deta) + \
                       (inverted_jacobi[2, 2]*shape_dpsi)
            # Изопараметрическая матрица градиентов
            b = zeros([6, 3*self.size])
            for j in range(0, self.size):
                b[0][j*3 + 0] = shape_dx[j]
                b[1][j*3 + 1] = shape_dy[j]
                b[2][j*3 + 2] = shape_dz[j]
                b[3][j*3 + 0] = shape_dy[j]
                b[3][j*3 + 1] = shape_dx[j]
                b[4][j*3 + 1] = shape_dz[j]
                b[4][j*3 + 2] = shape_dy[j]
                b[5][j*3 + 0] = shape_dz[j]
                b[5][j*3 + 2] = shape_dx[j]
            local_k += b.conj().transpose().dot(self._elastic_matrix_()).dot(b)*abs(jacobian)*self._w_[i]
        return local_k

    # Формирование локальных матриц масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        local_m = zeros((3*self.size, 3*self.size))
        # Интегрирование по кубу [-1; 1] x [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w_)):
            # Изопараметрические функции формы и их производные
            shape = self._shape_(i)
            shape_dxi = self._shape_dxi_(i)
            shape_deta = self._shape_deta_(i)
            shape_dpsi = self._shape_dpsi_(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1]), sum(shape_dxi*self.x[:, 2])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1]), sum(shape_deta*self.x[:, 2])],
                [sum(shape_dpsi*self.x[:, 0]), sum(shape_dpsi*self.x[:, 1]), sum(shape_dpsi*self.x[:, 2])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                [shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3], 0.0, 0.0,
                 shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7], 0.0, 0.0],
                [0.0, shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3], 0.0,
                 0.0, shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7], 0.0],
                [0.0, 0.0, shape[0], 0.0, 0.0, shape[1], 0.0, 0.0, shape[2], 0.0, 0.0, shape[2], 0.0, 0.0, shape[3],
                 0.0, 0.0, shape[4], 0.0, 0.0, shape[5], 0.0, 0.0, shape[6], 0.0, 0.0, shape[7]]
            ])
            local_m += c.conj().transpose().dot(c)*abs(jacobian)*self._w_[i]
        m = local_m*self.density
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Четырехугольный КЭ плиты
class TFE2D4P(TFE2D4):
    def __init__(self):
        super().__init__()

    def _extra_elastic_matrix_(self):
        return array([
            [1.0, 0.0],
            [0.0, 1.0],
        ])*self.e[0]/(2.0 + 2.0*self.m[0])

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b1 = zeros([3, 3*self.size])
            b2 = zeros([2, 3*self.size])
            for j in range(0, self.size):
                shape = self.a[j][0] + self.a[j][1]*self.x[i][0] + self.a[j][2]*self.x[i][1] + \
                        self.a[j][3]*self.x[i][0]*self.x[i][1]
                dx = self.a[j][1] + self.a[j][3]*self.x[i][1]
                dy = self.a[j][2] + self.a[j][3]*self.x[i][0]
                b1[0][3*j + 2] = -dx
                b1[1][3*j + 1] = dy
                b1[2][3*j + 1] = dx
                b1[2][3*j + 2] = -dy
                b2[0][3*j + 0] = dx
                b2[0][3*j + 2] = shape
                b2[1][3*j + 0] = dy
                b2[1][3*j + 1] = -shape
            e1 = b1.dot(u)
            s1 = self._elastic_matrix_().dot(e1)
            e2 = b2.dot(u)
            s2 = self._extra_elastic_matrix_().dot(e2)
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
    def _generate_stiffness_matrix_(self):
        local_k = zeros((3*self.size, 3*self.size))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w_)):
            # Изопараметрические функции формы и их производные
            shape = self._shape_(i)
            shape_dxi = self._shape_dxi_(i)
            shape_deta = self._shape_deta_(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta
            shape_dy = inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta
            # Изопараметрические матрицы градиентов
            b1 = zeros([3, 3*self.size])
            for j in range(0, self.size):
                b1[0][3*j + 2] = -shape_dx[j]
                b1[1][3*j + 1] = shape_dy[j]
                b1[2][3*j + 1] = shape_dx[j]
                b1[2][3*j + 2] = -shape_dy[j]
            b2 = zeros([2, 3*self.size])
            for j in range(0, self.size):
                b2[0][3*j + 0] = shape_dx[j]
                b2[0][3*j + 2] = shape[j]
                b2[1][3*j + 0] = shape_dy[j]
                b2[1][3*j + 1] = -shape[j]
            local_k += (b1.conj().transpose().dot(self._elastic_matrix_()).
                        dot(b1)*self.h**3/12.0 + b2.conj().transpose().dot(self._extra_elastic_matrix_()).
                        dot(b2)*self.h*5.0/6.0)*abs(jacobian)*self._w_[i]
        return local_k

    # Формирование локальных матриц жесткости, масс и демпфирования
    def _generate_mass_damping_matrix_(self):
        local_m = zeros((3*self.size, 3*self.size))
        # Интегрирование по прямоугольнику [-1; 1] x [-1; 1] (по формуле Гаусса)
        for i in range(len(self._w_)):
            # Изопараметрические функции формы и их производные
            shape = self._shape_(i)
            shape_dxi = self._shape_dxi_(i)
            shape_deta = self._shape_deta_(i)
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                [shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0, shape[3], 0, 0],
                [0, shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0, shape[3], 0],
                [0, 0, shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0, shape[3]],
            ])
            mi = array([
                [self.density*self.h, 0, 0],
                [0, self.density*self.h**3/12.0, 0],
                [0, 0, self.density*self.h**3/12.0],
            ])
            local_m += c.conj().transpose().dot(mi).dot(c)*abs(jacobian)*self._w_[i]
        m = local_m*self.density
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Треугольный КЭ пластины
class TFE2D3P(TFE2D3):
    def __init__(self):
        super().__init__()
        # Параметры квадратур Гаусса
        self._xi_ = [0, 1/2, 1/2]
        self._eta_ = [1/2, 0, 1/2]
        self._w_ = [1/6, 1/6, 1/6]

    def _extra_elastic_matrix_(self):
        # Вспомогательная матрица упругих свойст
        return array([
            [1.0, 0.0],
            [0.0, 1.0],
        ])*self.e[0]/(2.0 + 2.0*self.m[0])

    def calc(self, u):
        res = zeros((12, self.size))
        for i in range(0, self.size):
            # Матрица градиентов
            b1 = zeros([3, 3*self.size])
            b2 = zeros([2, 3*self.size])
            for j in range(0, self.size):
                shape = self.a[j][0] + self.a[j][1]*self.x[i][0] + self.a[j][2]*self.x[i][1]
                dx = self.a[j][1]
                dy = self.a[j][2]
                b1[0][3*j + 2] = -dx
                b1[1][3*j + 1] = dy
                b1[2][3*j + 1] = dx
                b1[2][3*j + 2] = -dy

                b2[0][3*j + 0] = dx
                b2[0][3*j + 2] = shape
                b2[1][3*j + 0] = dy
                b2[1][3*j + 1] = -shape
            e1 = b1.dot(u)
            s1 = self._elastic_matrix_().dot(e1)
            e2 = b2.dot(u)
            s2 = self._extra_elastic_matrix_().dot(e2)
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

    def _generate_stiffness_matrix_(self):
        local_k = zeros((3*self.size, 3*self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[0,1] (по формуле Гаусса)
        for i in range(len(self._w_)):
            # Изопараметрические функции формы и их производные
            shape = array([1.0 - self._xi_[i] - self._eta_[i], self._xi_[i], self._eta_[i]])
            shape_dxi = array([-1.0, 1.0, 0.0])
            shape_deta = array([-1.0, 0.0, 1.0])
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            inverted_jacobi = inv(jacobi)
            shape_dx = inverted_jacobi[0, 0]*shape_dxi + inverted_jacobi[0, 1]*shape_deta
            shape_dy = inverted_jacobi[1, 0]*shape_dxi + inverted_jacobi[1, 1]*shape_deta
            # Матрицы градиентов
            b1 = zeros([3, 3*self.size])
            for j in range(0, self.size):
                b1[0][3*j + 2] = -shape_dx[j]
                b1[1][3*j + 1] = shape_dy[j]
                b1[2][3*j + 1] = shape_dx[j]
                b1[2][3*j + 2] = -shape_dy[j]
            b2 = zeros([2, 3*self.size])
            for j in range(0, self.size):
                b2[0][3*j + 0] = shape_dx[j]
                b2[0][3*j + 2] = shape[j]
                b2[1][3*j + 0] = shape_dy[j]
                b2[1][3*j + 1] = -shape[j]
            # Вычисление компонент локальной матрицы жесткости
            local_k += (b1.conj().transpose().dot(self._elastic_matrix_()).
                        dot(b1)*self.h**3/12.0 + b2.conj().transpose().
                        dot(self._extra_elastic_matrix_()).dot(b2)*self.h*5.0/6.0)*abs(jacobian)*self._w_[i]
        return local_k

    def _generate_mass_damping_matrix_(self):
        local_m = zeros((3*self.size, 3*self.size))
        # Интегрирование по треугольнику [0,0]-[1,0]-[1,1] (по формуле Гаусса)
        for i in range(len(self._w_)):
            # Изопараметрические функции формы и их производные
            shape = array([1.0 - self._xi_[i] - self._eta_[i], self._xi_[i], self._eta_[i]])
            shape_dxi = array([-1.0, 1.0, 0.0])
            shape_deta = array([0.0, -1.0, 1.0])
            # Матрица Якоби
            jacobi = array([
                [sum(shape_dxi*self.x[:, 0]), sum(shape_dxi*self.x[:, 1])],
                [sum(shape_deta*self.x[:, 0]), sum(shape_deta*self.x[:, 1])]
            ])
            # Якобиан
            jacobian = det(jacobi)
            # Вспомогательная матрица для построения матриц масс и демпфирования
            c = array([
                [shape[0], 0, 0, shape[1], 0, 0, shape[2], 0, 0],
                [0, shape[0], 0, 0, shape[1], 0, 0, shape[2], 0],
                [0, 0, shape[0], 0, 0, shape[1], 0, 0, shape[2]]
            ])
            mi = array([
                [self.h, 0, 0],
                [0, self.h**3/12.0, 0],
                [0, 0, self.h**3/12.0],
            ])
            local_m += c.conj().transpose().dot(mi).dot(c)*abs(jacobian)*self._w_[i]
        m = local_m*self.density
        c = m*self.damping[0] + self.K*self.damping[1]
        return m, c


# Треугольный КЭ оболочки
class TFE2D3S(TFE2D3P):
    def __init__(self):
        super().__init__()
        self.size = 3
        self.T = zeros((3, 3))  # Матрица преобразования глобальных координат в локальные
        self.global_x = zeros((3, 3))

    def _create_(self):
        self.T = self._create_transform_matrix_()
        self.global_x = self.x
        self.x = array([self.T.dot(self.x[0, :]), self.T.dot(self.x[1, :]), self.T.dot(self.x[2, :])])
        self.x -= self.x[0, :]
        TFE2D3._create_(self)
#        self.x = x

    def _prepare_transform_matrix_(self):
        m = zeros((6*self.size, 6*self.size))
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
            membrane_u.append(lu[6*i])      # u
            membrane_u.append(lu[6*i + 1])  # w
        membrane_res = TFE2D3.calc(self, membrane_u)
        # Вычисляем деформации и напряжения пластины
        plate_u = []
        for i in range(0, self.size):
            plate_u.append(lu[6*i + 2])     # w
            plate_u.append(lu[6*i + 3])     # Tx
            plate_u.append(lu[6*i + 4])     # Ty
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

    def _generate_stiffness_matrix_(self):
        local_k = identity(6*self.size)
        # Создание матрицы преобразования
        m = self._prepare_transform_matrix_()
        # Локальная матрицыа жесткости плоского КЭ
        k1 = TFE2D3._generate_stiffness_matrix_(self)
        # ... КЭ пластины
        k2 = TFE2D3P._generate_stiffness_matrix_(self)
        global_freedom = 6
        local_freedom1 = 2
        for i in range(0, len(k1)):
            p = (i//local_freedom1)*global_freedom + i % local_freedom1
            for j in range(0, len(k1)):
                q = (j//local_freedom1)*global_freedom + j % local_freedom1
                local_k[p][q] = k1[i][j]
        local_freedom2 = 3
        for i in range(0, len(k2)):
            p = (i//local_freedom2)*global_freedom + i % local_freedom2 + local_freedom1
            for j in range(0, len(k2)):
                q = (j//local_freedom2)*global_freedom + j % local_freedom2 + local_freedom1
                local_k[p][q] = k2[i][j]
        return m.conj().transpose().dot(local_k).dot(m)


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
            membrane_u.append(lu[6*i])      # u
            membrane_u.append(lu[6*i + 1])  # w
        membrane_res = TFE2D4.calc(self, membrane_u)
        # Вычисляем деформации и напряжения пластины
        plate_u = []
        for i in range(0, self.size):
            plate_u.append(lu[6*i + 2])     # w
            plate_u.append(lu[6*i + 3])     # Tx
            plate_u.append(lu[6*i + 4])     # Ty
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
        m = zeros((6*self.size, 6*self.size))
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
    def _generate_stiffness_matrix_(self):
        local_k = identity(6*self.size)
        # Подготовка матрицы преобразования
        m = self._prepare_transform_matrix_()
        # Локальная матрицыа жесткости плоского КЭ
        k1 = TFE2D4._generate_stiffness_matrix_(self)

        # ... КЭ пластины
        k2 = TFE2D4P._generate_stiffness_matrix_(self)
        global_freedom = 6
        local_freedom1 = 2
        for i in range(0, len(k1)):
            p = (i//local_freedom1)*global_freedom + i % local_freedom1
            for j in range(0, len(k1)):
                q = (j//local_freedom1)*global_freedom + j % local_freedom1
                local_k[p][q] = k1[i][j]
        local_freedom2 = 3
        for i in range(0, len(k2)):
            p = (i//local_freedom2)*global_freedom + i % local_freedom2 + local_freedom1
            for j in range(0, len(k2)):
                q = (j//local_freedom2)*global_freedom + j % local_freedom2 + local_freedom1
                local_k[p][q] = k2[i][j]

        # import sys
        # print('\n******************************************')
        # for i in range(0, len(local_k)):
        #     for j in range(0, len(local_k)):
        #         sys.stdout.write('%+E\t' % local_k[i][j])
        #     sys.stdout.write('\n')
        # print('******************************************')
        return m.conj().transpose().dot(local_k).dot(m)

    def _create_(self):
        self.T = self._create_transform_matrix_()
        self.global_x = self.x

        x0 = self.T.dot(self.x[0, :].conj().transpose())
        x1 = self.T.dot(self.x[1, :].conj().transpose())
        x2 = self.T.dot(self.x[2, :].conj().transpose())
        x3 = self.T.dot(self.x[3, :].conj().transpose())

        self.x = array([x0 - x0, x1 - x0, x2 - x0, x3 - x0])
#        self.x = array([x0, x1, x2, x3])
        TFE2D4._create_(self)
#        self.x = x
