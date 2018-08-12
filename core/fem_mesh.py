#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Конечно-элементная модель объекта расчета
###################################################################

import math
from core.fem_error import TFEMException

# Типы конечных элементов
FEType = [
    'fe_1d_2',
    'fe_2d_3',
    'fe_2d_3_p',
    'fe_2d_3_s',
    'fe_2d_4',
    'fe_2d_4_p',
    'fe_3d_4',
    'fe_3d_8'
]


class TMesh:
    def __init__(self):
        self.mesh_file = ''     # Имя файла с данными о геометрической модели
        self.fe_type = ''       # Тип КЭ
        self.x = []             # Координаты узлов
        self.be = []            # Связи граничных элементов
        self.fe = []            # Связи в КЭ
        self.freedom = 0        # Кол-во степеней свободы
        self.dimension = 0      # Размерность КЭ

    @staticmethod
    def get_fe_data(t):
        if t == 3:
            return 'fe_2d_3', 2, 3, 2, 2
        elif t == 4:
            return 'fe_3d_4', 3, 4, 3, 3
        elif t == 8:
            return 'fe_3d_8', 4, 8, 3, 3
        elif t == 24:
            return 'fe_2d_4', 2, 4, 2, 2
        elif t == 34:
            return 'fe_1d_2', 0, 2, 1, 1
        elif t == 123:
            return 'fe_2d_3_p', 0, 3, 3, 2
        elif t == 124:
            return 'fe_2d_4_p', 0, 4, 3, 2
        elif t == 223:
            return 'fe_2d_3_s', 0, 3, 6, 3
        elif t == 224:
            return 'fe_2d_4_s', 0, 4, 6, 3
        else:
            raise TFEMException('unknown_fe_err')

    def load(self, name):
        try:
            self.mesh_file = name
            file = open(self.mesh_file)
            lines = file.readlines()
            file.close()
        except IOError:
            raise TFEMException('read_file_err')
        self.fe_type, size_surface, size_fe, self.freedom, self.dimension = self.get_fe_data(int(lines[0]))
        # Кол-во узлов
        n = int(lines[1])
        # Считываем узлы
        index = 2
        for i in range(0, n):
            tx = list()
            tx.append(float(lines[2 + i].split()[0]))
            if self.dimension > 1:
                tx.append(float(lines[2 + i].split()[1]))
            if self.dimension > 2:
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
        for i in range(0, n):
            row = []
            for j in range(0, size_surface):
                row.append(int(lines[index].split()[j]))
            self.be.append(row)
            index += 1
        if self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_3_s' or self.fe_type == 'fe_2d_4_p' or \
                self.fe_type == 'fe_2d_4_s':
            self.be = self.fe

    def fe_name(self):
        if self.fe_type == 'fe_1d_2':
            return 'one-dimensional linear element (2 nodes)'
        elif self.fe_type == 'fe_2d_3':
            return 'linear triangular element (3 nodes)'
        elif self.fe_type == 'fe_2d_4':
            return 'quadrilateral element (4 nodes)'
        elif self.fe_type == 'fe_2d_3_p':
            return 'plate element (3 nodes)'
        elif self.fe_type == 'fe_2d_3_s':
            return 'shell element (3 nodes)'
        elif self.fe_type == 'fe_2d_4_p':
            return 'plate element (4 nodes)'
        elif self.fe_type == 'fe_2d_4_s':
            return 'shell element (4 nodes)'
        elif self.fe_type == 'fe_3d_4':
            return 'linear tetrahedron (4 nodes)'
        elif self.fe_type == 'fe_3d_8':
            return 'cube element (8 nodes)'

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

    def get_be_coord(self, index):
        x = []
        for i in range(0, len(self.be[index])):
            a = []
            for j in range(0, self.dimension):
                a.append(self.x[self.be[index][i]][j])
            x.append(a)
        return x

    # Вычисление длины (площади) заданного граничного элемента
    def square(self, index):
        x = self.get_be_coord(index)
        s = 0
        if self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_4':  # ГЭ - отрезок
            s = math.sqrt((x[0][0] - x[1][0])**2 + (x[0][1] - x[1][1])**2)
        elif self.fe_type == 'fe_2d_3_p':  # ГЭ - "плоский" треугольник
            a = math.sqrt((x[0][0] - x[1][0])**2 + (x[0][1] - x[1][1])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2)
            c = math.sqrt((x[2][0] - x[1][0])**2 + (x[2][1] - x[1][1])**2)
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif self.fe_type == 'fe_2d_3_s' or self.fe_type == 'fe_3d_4':  # ГЭ - "пространственный" треугольник
            a = math.sqrt((x[0][0] - x[1][0])**2 + (x[0][1] - x[1][1])**2 + (x[0][2] - x[1][2])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2 + (x[0][2] - x[2][2])**2)
            c = math.sqrt((x[2][0] - x[1][0])**2 + (x[2][1] - x[1][1])**2 + (x[2][2] - x[1][2])**2)
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif self.fe_type == 'fe_2d_4_p':  # ГЭ - "плоский" четырехугольник
            a = math.sqrt((x[0][0] - x[1][0])**2 + (x[0][1] - x[1][1])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2)
            c = math.sqrt((x[2][0] - x[1][0])**2 + (x[2][1] - x[1][1])**2)
            p = 0.5*(a + b + c)
            s = math.sqrt(p * (p - a) * (p - b) * (p - c))
            a = math.sqrt((x[0][0] - x[3][0])**2 + (x[0][1] - x[3][1])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2)
            c = math.sqrt((x[2][0] - x[3][0])**2 + (x[2][1] - x[3][1])**2)
            p = 0.5*(a + b + c)
            s += math.sqrt(p * (p - a) * (p - b) * (p - c))
        elif self.fe_type == 'fe_3d_8' or self.fe_type == 'fe_2d_4_s':   # ГЭ - четырехугольник
            a = math.sqrt((x[0][0] - x[1][0])**2 + (x[0][1] - x[1][1])**2 + (x[0][2] - x[1][2])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2 + (x[0][2] - x[2][2])**2)
            c = math.sqrt((x[2][0] - x[1][0])**2 + (x[2][1] - x[1][1])**2 + (x[2][2] - x[1][2])**2)
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))
            a = math.sqrt((x[0][0] - x[3][0])**2 + (x[0][1] - x[3][1])**2 + (x[0][2] - x[3][2])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2 + (x[0][2] - x[2][2])**2)
            c = math.sqrt((x[2][0] - x[3][0])**2 + (x[2][1] - x[3][1])**2 + (x[2][2] - x[3][2])**2)
            p = 0.5*(a + b + c)
            s += math.sqrt(p*(p - a)*(p - b)*(p - c))
        return s

    # Вычисление объема (длины, площади) заданного конечного элемента
    def volume(self, index):
        v = 0
        x = self.get_fe_coord(index)
        if self.fe_type == 'fe_1d_2':
            v = math.fabs(x[0][0] - x[1][0])
        elif self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_3_s':
            a = math.sqrt((x[0][0] - x[1][0])**2 + (x[0][1] - x[1][1])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2)
            c = math.sqrt((x[2][0] - x[1][0])**2 + (x[2][1] - x[1][1])**2)
            p = 0.5*(a + b + c)
            v = math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif self.fe_type == 'fe_2d_4' or self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_2d_4_s':
            a = math.sqrt((x[0][0] - x[1][0])**2 + (x[0][1] - x[1][1])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2)
            c = math.sqrt((x[2][0] - x[1][0])**2 + (x[2][1] - x[1][1])**2)
            p = 0.5*(a + b + c)
            v = math.sqrt(p*(p - a)*(p - b)*(p - c))
            a = math.sqrt((x[0][0] - x[3][0])**2 + (x[0][1] - x[3][1])**2)
            b = math.sqrt((x[0][0] - x[2][0])**2 + (x[0][1] - x[2][1])**2)
            c = math.sqrt((x[2][0] - x[3][0])**2 + (x[2][1] - x[3][1])**2)
            p = 0.5*(a + b + c)
            v += math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif self.fe_type == 'fe_3d_4':
            a = (x[1][0] - x[0][0]) * (x[2][1] - x[0][1]) * (x[3][2] - x[0][2]) + (x[3][0] - x[0][0]) * \
                (x[1][1] - x[0][1]) * (x[2][2] - x[0][2]) + (x[2][0] - x[0][0])*(x[3][1] - x[0][1]) * \
                (x[1][2] - x[0][2])
            b = (x[3][0] - x[0][0]) * (x[2][1] - x[0][1]) * (x[1][2] - x[0][2]) + (x[2][0] - x[0][0]) * \
                (x[1][1] - x[0][1]) * (x[3][2] - x[0][2]) + (x[1][0] - x[0][0])*(x[3][1] - x[0][1]) * \
                (x[2][2] - x[0][2])
            v = math.fabs(a - b)/6.0
        elif self.fe_type == 'fe_3d_8':
            ref = [[0, 1, 4, 7], [4, 1, 5, 7], [1, 2, 6, 7], [1, 5, 6, 7], [1, 2, 3, 7], [0, 3, 1, 7]]
            for i in range(0, 6):
                a = (x[ref[i][1]][0] - x[ref[i][0]][0]) * (x[ref[i][2]][1] - x[ref[i][0]][1]) * \
                    (x[ref[i][3]][2] - x[ref[i][0]][2]) + (x[ref[i][3]][0] - x[ref[i][0]][0]) * \
                    (x[ref[i][1]][1] - x[ref[i][0]][1])*(x[ref[i][2]][2] - x[ref[i][0]][2]) + \
                    (x[ref[i][2]][0] - x[ref[i][0]][0]) * (x[ref[i][3]][1] - x[ref[i][0]][1]) * \
                    (x[ref[i][1]][2] - x[ref[i][0]][2])
                b = (x[ref[i][3]][0] - x[ref[i][0]][0]) * (x[ref[i][2]][1] - x[ref[i][0]][1]) * \
                    (x[ref[i][1]][2] - x[ref[i][0]][2]) + (x[ref[i][2]][0] - x[ref[i][0]][0]) * \
                    (x[ref[i][1]][1] - x[ref[i][0]][1]) * (x[ref[i][3]][2] - x[ref[i][0]][2]) + \
                    (x[ref[i][1]][0] - x[ref[i][0]][0]) * (x[ref[i][3]][1] - x[ref[i][0]][1]) * \
                    (x[ref[i][2]][2] - x[ref[i][0]][2])
                v += math.fabs(a - b)/6.0
        return v

    def is_plate(self):
        return True if self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_4_p' else False

    def is_shell(self):
        return True if self.fe_type == 'fe_2d_3_s' or self.fe_type == 'fe_2d_4_s' else False
