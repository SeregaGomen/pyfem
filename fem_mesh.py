#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Конечно-элементная модель объекта расчета
###################################################################

from fem_error import TFEMException

# Типы конечных элементов
FEType = [
    'fe_1d_2',
    'fe_2d_3',
    'fe_2d_4',
    'fe_3d_4',
    'fe_3d_8'
]


class TMesh:
    def __init__(self):
        self.mesh_file = ''     # Имя файла с данными о геометрической модели
        self.fe_type = ''       # Тип КЭ
        self.x = []             # Координаты узлов
        self.y = []
        self.z = []
        self.surface = []       # Связи граничных элементов
        self.fe = []            # Связи в КЭ
        self.freedom = 0        # Кол-во степеней свободы

    @staticmethod
    def get_fe_data(t):
        if t == 3:
            return 'fe_2d_3', 2, 3, 2
        elif t == 4:
            return 'fe_3d_4', 3, 4, 3
        elif t == 8:
            return 'fe_3d_8', 4, 8, 3
        elif t == 24:
            return 'fe_2d_4', 2, 4, 2
        elif t == 34:
            return 'fe_1d_2', 0, 2, 1
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
        self.fe_type, size_surface, size_fe, self.freedom = self.get_fe_data(int(lines[0]))
        # Кол-во узлов
        n = int(lines[1])
        # Считываем узлы
        index = 2
        for i in range(0, n):
            self.x.append(float(lines[2 + i].split()[0]))
            if self.freedom > 1:
                self.y.append(float(lines[2 + i].split()[1]))
            if self.freedom > 2:
                self.z.append(float(lines[2 + i].split()[2]))
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
            self.surface.append(row)
            index += 1

    def fe_name(self, t):
        if self.fe_type == 'fe_1d_2':
            return 'one-dimensional linear element (2 nodes)'
        elif self.fe_type == 'fe_2d_3':
            return 'linear triangular element (3 nodes)'
        elif self.fe_type == 'fe_2d_4':
            return 'quadrilateral element (4 nodes)'
        elif self.fe_type == 'fe_2d_6':
            return 'quadratic triangular element (6 nodes)'
        elif self.fe_type == 'fe_3d_4':
            return 'linear tetrahedron (4 nodes)'
        elif self.fe_type == 'fe_3d_8':
            return 'cube element (8 nodes)'
        elif self.fe_type == 'fe_3d_10':
            return 'quadratic tetrahedron (10 nodes)'

    def get_coord(self, i):
        return self.x[i], self.y[i] if (len(self.y)) else 0, self.z[i] if (len(self.z)) else 0
