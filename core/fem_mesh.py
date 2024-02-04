#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Конечно-элементная модель объекта расчета
###################################################################

from abc import abstractmethod
from math import sqrt
from numpy import array
from numpy.linalg import det
from core.fem_error import TException


# Типы конечных элементов
FEType = [
    'fe_1d_2',
    'fe_2d_3',
    'fe_2d_4',
    'fe_2d_6',
    'fe_2d_3_p',
    'fe_2d_4_p',
    'fe_2d_6_p',
    'fe_3d_3_s',
    'fe_3d_4_s',
    'fe_3d_6_s',
    'fe_3d_4',
    'fe_3d_8',
    'fe_3d_10'
]


# Абстрактный класс, инкапсулирующий понятие дискретной модели объекта
class TMesh:
    def __init__(self):
        self.mesh_file = ''         # Имя файла с данными о геометрической модели
        self.fe_type = ''           # Тип КЭ
        self.x = []                 # Координаты узлов
        self.be = []                # Связи граничных элементов
        self.fe = []                # Связи в КЭ
        self.freedom = 0            # Кол-во степеней свободы
        self.dimension = 0          # Размерность КЭ

    @staticmethod
    def get_fe_data(t):
        if t == '3' or t == 'fe2d3':
            return 'fe_2d_3', 2, 3, 2, 2
        elif t == '4' or t == 'fe2d34':
            return 'fe_3d_4', 3, 4, 3, 3
        elif t == '6' or t == 'fe2d6':
            return 'fe_2d_6', 3, 6, 2, 2
        elif t == '8' or t == 'fe3d8':
            return 'fe_3d_8', 4, 8, 3, 3
        elif t == '10' or t == 'fe3d10':
            return 'fe_3d_10', 6, 10, 3, 3
        elif t == '24' or t == 'fe2d4':
            return 'fe_2d_4', 2, 4, 2, 2
        elif t == '34' or t == 'fe1d2':
            return 'fe_1d_2', 1, 2, 1, 1
        elif t == '123' or t == 'fe2d3p':
            return 'fe_2d_3_p', 0, 3, 3, 2
        elif t == '124' or t == 'fe2d4p':
            return 'fe_2d_4_p', 0, 4, 3, 2
        elif t == '125' or t == 'fe2d6p':
            return 'fe_2d_6_p', 0, 6, 3, 2
        elif t == '223' or t == 'fe3d3s':
            return 'fe_3d_3_s', 0, 3, 6, 3
        elif t == '224' or t == 'fe3d4s':
            return 'fe_3d_4_s', 0, 4, 6, 3
        elif t == '225' or t == 'fe3d6s':
            return 'fe_3d_6_s', 0, 6, 6, 3
        else:
            raise TException('unknown_fe_err')

    @abstractmethod
    def load(self, name):
        raise NotImplementedError('Method TMesh.load() is pure virtual')

    def fe_name(self):
        if self.fe_type == 'fe_1d_2':
            return 'one-dimensional linear element (2 nodes)'
        elif self.fe_type == 'fe_2d_3':
            return 'linear triangular element (3 nodes)'
        elif self.fe_type == 'fe_2d_4':
            return 'quadrilateral element (4 nodes)'
        elif self.fe_type == 'fe_2d_6':
            return 'quadratic triangular element (6 nodes)'
        elif self.fe_type == 'fe_2d_3_p':
            return 'triangular plate element (3 nodes)'
        elif self.fe_type == 'fe_2d_6_p':
            return 'triangular plate element (6 nodes)'
        elif self.fe_type == 'fe_3d_3_s':
            return 'triangular shell element (3 nodes)'
        elif self.fe_type == 'fe_3d_6_s':
            return 'triangular shell element (6 nodes)'
        elif self.fe_type == 'fe_2d_4_p':
            return 'plate element (4 nodes)'
        elif self.fe_type == 'fe_3d_4_s':
            return 'shell element (4 nodes)'
        elif self.fe_type == 'fe_3d_4':
            return 'linear tetrahedron (4 nodes)'
        elif self.fe_type == 'fe_3d_8':
            return 'cube element (8 nodes)'
        elif self.fe_type == 'fe_3d_10':
            return 'quadratic tetrahedron (10 nodes)'

    def get_coord(self, index):
        return self.x[index]

    def get_fe_coord(self, index):
        x = []
        for i in range(0, len(self.fe[index])):
            a = []
            for j in range(0, self.dimension):
                a.append(self.x[self.fe[index][i]][j])
            x.append(a)
        return x

    def get_fe_center(self, index):
        x = []
        for j in range(self.dimension):
            a = 0
            for i in range(len(self.fe[index])):
                a += self.x[self.fe[index][i]][j]
            x.append(a / self.dimension)
        return x

    def get_be_coord(self, index):
        x = []
        for i in range(0, len(self.be[index])):
            a = []
            for j in range(0, self.dimension):
                a.append(self.x[self.be[index][i]][j])
            x.append(a)
        return x

    def get_be_center(self, index):
        x = []
        for j in range(self.dimension):
            a = 0
            for i in range(len(self.be[index])):
                a += self.x[self.be[index][i]][j]
            x.append(a / self.dimension)
        return x

    def is_plate(self):
        return True if self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_2d_4_p' \
            else False

    def is_shell(self):
        return True if self.fe_type == 'fe_3d_3_s' or self.fe_type == 'fe_3d_4_s' or self.fe_type == 'fe_3d_6_s' \
            else False

    def is_1d(self):
        return True if self.fe_type == 'fe_1d_2' else False

    def is_2d(self):
        return True if self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_4' or self.fe_type == 'fe_2d_6' else False

    def is_3d(self):
        return True if self.fe_type == 'fe_3d_4' or self.fe_type == 'fe_3d_8' or self.fe_type == 'fe_3d_10' else False

    def base_be_size(self):
        if self.fe_type == 'fe_2d_6':
            return 2
        elif self.fe_type == 'fe_3d_10' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_3d_6_s':
            return 3
        return len(self.be[0])

    def base_fe_size(self):
        if self.fe_type == 'fe_2d_6' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_3d_6_s':
            return 3
        elif self.fe_type == 'fe_3d_10':
            return 4
        return len(self.fe[0])

    # Вычисление длины (площади) заданного граничного элемента
    def square(self, index):
        x = self.get_be_coord(index)
        if self.dimension == 2:     # Добавление z = 0 для плоского случая
            for i in range(0, len(x)):
                x[i].append(0)
        s = 0
        if self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_4' or self.fe_type == 'fe_2d_6':  # ГЭ - отрезок
            s = sqrt((x[0][0] - x[1][0]) ** 2 + (x[0][1] - x[1][1]) ** 2)
        elif self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_3d_3_s' or self.fe_type == 'fe_3d_4' or \
                self.fe_type == 'fe_3d_10' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_3d_6_s':
            # ГЭ - треугольник
            a = sqrt((x[0][0] - x[1][0]) ** 2 + (x[0][1] - x[1][1]) ** 2 + (x[0][2] - x[1][2]) ** 2)
            b = sqrt((x[0][0] - x[2][0]) ** 2 + (x[0][1] - x[2][1]) ** 2 + (x[0][2] - x[2][2]) ** 2)
            c = sqrt((x[2][0] - x[1][0]) ** 2 + (x[2][1] - x[1][1]) ** 2 + (x[2][2] - x[1][2]) ** 2)
            p = 0.5 * (a + b + c)
            s = sqrt(p * (p - a) * (p - b) * (p - c))
        elif self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_3d_8' or self.fe_type == 'fe_3d_4_s':
            # ГЭ - четырехугольник
            a = sqrt((x[0][0] - x[1][0]) ** 2 + (x[0][1] - x[1][1]) ** 2 + (x[0][2] - x[1][2]) ** 2)
            b = sqrt((x[0][0] - x[2][0]) ** 2 + (x[0][1] - x[2][1]) ** 2 + (x[0][2] - x[2][2]) ** 2)
            c = sqrt((x[2][0] - x[1][0]) ** 2 + (x[2][1] - x[1][1]) ** 2 + (x[2][2] - x[1][2]) ** 2)
            p = 0.5 * (a + b + c)
            s = sqrt(p * (p - a) * (p - b) * (p - c))
            a = sqrt((x[0][0] - x[3][0]) ** 2 + (x[0][1] - x[3][1]) ** 2 + (x[0][2] - x[3][2]) ** 2)
            b = sqrt((x[0][0] - x[2][0]) ** 2 + (x[0][1] - x[2][1]) ** 2 + (x[0][2] - x[2][2]) ** 2)
            c = sqrt((x[2][0] - x[3][0]) ** 2 + (x[2][1] - x[3][1]) ** 2 + (x[2][2] - x[3][2]) ** 2)
            p = 0.5 * (a + b + c)
            s += sqrt(p * (p - a) * (p - b) * (p - c))
        return s

    # Вычисление объема (длины, площади) заданного конечного элемента
    def volume(self, index):
        v = 0
        x = self.get_fe_coord(index)
        if self.fe_type == 'fe_1d_2':
            v = abs(x[0][0] - x[1][0])
        elif self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_6' or self.fe_type == 'fe_2d_3_p' or \
                self.fe_type == 'fe_3d_3_s' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_3d_6_s':
            a = sqrt((x[0][0] - x[1][0]) ** 2 + (x[0][1] - x[1][1]) ** 2)
            b = sqrt((x[0][0] - x[2][0]) ** 2 + (x[0][1] - x[2][1]) ** 2)
            c = sqrt((x[2][0] - x[1][0]) ** 2 + (x[2][1] - x[1][1]) ** 2)
            p = 0.5 * (a + b + c)
            v = sqrt(p * (p - a) * (p - b) * (p - c))
        elif self.fe_type == 'fe_2d_4' or self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_3d_4_s':
            a = sqrt((x[0][0] - x[1][0]) ** 2 + (x[0][1] - x[1][1]) ** 2)
            b = sqrt((x[0][0] - x[2][0]) ** 2 + (x[0][1] - x[2][1]) ** 2)
            c = sqrt((x[2][0] - x[1][0]) ** 2 + (x[2][1] - x[1][1]) ** 2)
            p = 0.5 * (a + b + c)
            v = sqrt(p * (p - a) * (p - b) * (p - c))
            a = sqrt((x[0][0] - x[3][0]) ** 2 + (x[0][1] - x[3][1]) ** 2)
            b = sqrt((x[0][0] - x[2][0]) ** 2 + (x[0][1] - x[2][1]) ** 2)
            c = sqrt((x[2][0] - x[3][0]) ** 2 + (x[2][1] - x[3][1]) ** 2)
            p = 0.5 * (a + b + c)
            v += sqrt(p * (p - a) * (p - b) * (p - c))
        elif self.fe_type == 'fe_3d_4' or self.fe_type == 'fe_3d_10':
            a = (x[1][0] - x[0][0]) * (x[2][1] - x[0][1]) * (x[3][2] - x[0][2]) + (x[3][0] - x[0][0]) * \
                (x[1][1] - x[0][1]) * (x[2][2] - x[0][2]) + (x[2][0] - x[0][0])*(x[3][1] - x[0][1]) * \
                (x[1][2] - x[0][2])
            b = (x[3][0] - x[0][0]) * (x[2][1] - x[0][1]) * (x[1][2] - x[0][2]) + (x[2][0] - x[0][0]) * \
                (x[1][1] - x[0][1]) * (x[3][2] - x[0][2]) + (x[1][0] - x[0][0])*(x[3][1] - x[0][1]) * \
                (x[2][2] - x[0][2])
            v = abs(a - b)/6.0
        elif self.fe_type == 'fe_3d_8':
            ref = [[0, 1, 4, 7], [4, 1, 5, 7], [1, 2, 6, 7], [1, 5, 6, 7], [1, 2, 3, 7], [0, 3, 1, 7]]
            for i in range(0, 6):
                m = array([
                    [x[ref[i][1]][0] - x[ref[i][0]][0], x[ref[i][1]][1] - x[ref[i][0]][1],
                     x[ref[i][1]][2] - x[ref[i][0]][2]],
                    [x[ref[i][2]][0] - x[ref[i][0]][0], x[ref[i][2]][1] - x[ref[i][0]][1],
                     x[ref[i][2]][2] - x[ref[i][0]][2]],
                    [x[ref[i][3]][0] - x[ref[i][0]][0], x[ref[i][3]][1] - x[ref[i][0]][1],
                     x[ref[i][3]][2] - x[ref[i][0]][2]],
                ])
                v += abs(det(m)) / 6
        return v

    # Построение вектора нормали к заданному ГЭ
    def be_normal(self, index):
        v = [0, 0, 0]
        # Построение нормали зависит от типа КЭ
        if self.is_1d():
            # Одномерная задача
            v[0] = 1
        elif self.is_2d():
            # Плоская задача
            v[0] = self.x[self.be[index][0]][1] - self.x[self.be[index][1]][1]
            v[1] = self.x[self.be[index][1]][0] - self.x[self.be[index][0]][0]
        elif self.is_plate():
            # Пластины
            v[0] = v[1] = 0
        else:
            # Оболочка или трехмерная задача
            v[0] = (self.x[self.be[index][1]][1] - self.x[self.be[index][0]][1]) * (self.x[self.be[index][2]][2] -
                    self.x[self.be[index][0]][2]) - (self.x[self.be[index][2]][1] - self.x[self.be[index][0]][1]) * \
                   (self.x[self.be[index][1]][2] - self.x[self.be[index][0]][2])
            v[1] = (self.x[self.be[index][2]][0] - self.x[self.be[index][0]][0]) * (self.x[self.be[index][1]][2] -
                    self.x[self.be[index][0]][2]) - (self.x[self.be[index][1]][0] - self.x[self.be[index][0]][0]) * \
                   (self.x[self.be[index][2]][2] - self.x[self.be[index][0]][2])
            v[2] = (self.x[self.be[index][1]][0] - self.x[self.be[index][0]][0]) * (self.x[self.be[index][2]][1] -
                    self.x[self.be[index][0]][1]) - (self.x[self.be[index][2]][0] - self.x[self.be[index][0]][0]) * \
                   (self.x[self.be[index][1]][1] - self.x[self.be[index][0]][1])
        # Нормализируем вектор
        len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5
        v[0] /= len
        v[1] /= len
        v[2] /= len
        return v


# Дискретная модель в формате TRPA
class TMeshTRPA(TMesh):
    def __init__(self):
        super().__init__()

    def load(self, name):
        try:
            self.mesh_file = name
            file = open(self.mesh_file)
            lines = file.readlines()
            file.close()
        except IOError:
            raise TException('read_file_err')
        self.fe_type, size_surface, size_fe, self.freedom, self.dimension = self.get_fe_data(lines[0].strip())
        # Кол-во узлов
        n = int(lines[1])
        # Считываем узлы
        index = 2
        for i in range(0, n):
            tx = list()
            tx.append(float(lines[2 + i].split()[0]))
            if self.dimension > 1:
                tx.append(float(lines[2 + i].split()[1]))
            if self.dimension > 2 and not self.is_plate():
                tx.append(float(lines[2 + i].split()[2]))
            self.x.append(tx)
            index += 1
        # Кол-во КЭ
        n = int(lines[index])
        index += 1
        # Считываем КЭ
        for i in range(0, n):
            row = []
            for j in range(0, size_fe):
                row.append(int(lines[index].split()[j]))
            self.fe.append(row)
            index += 1
        # Кол-во ГЭ
        n = int(lines[index])
        index += 1
        # Считываем ГЭ
        for i in range(n):
            row = []
            for j in range(size_surface):
                row.append(int(lines[index].split()[j]))
            self.be.append(row)
            index += 1
        if self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_3d_3_s' or \
                self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_3d_6_s' or self.fe_type == 'fe_3d_4_s':
            self.be = self.fe
