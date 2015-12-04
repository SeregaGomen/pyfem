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
    'fe_2d_6',
    'fe_3d_4',
    'fe_3d_8',
    'fe_3d_10'
]


class TMesh:
    def __init__(self):
        self.fe_type = ''   # Тип КЭ
        self.x = []         # Координаты узлов
        self.y = []
        self.z = []
        self.surface = []   # Связи граничных элементов
        self.fe = []        # Связи в КЭ
        self.freedom = 0    # Кол-во степеней свободы

    def load(self, name):
        try:
            file = open(name)
            lines = file.readlines()
            file.close()
        except IOError:
            raise TFEMException('read_file_err')
        t = int(lines[0])
        if t == 3:
            self.fe_type = 'fe_2d_3'
            size_surface = 2
            size_fe = 3
            self.freedom = 2
        elif t == 4:
            self.fe_type = 'fe_3d_4'
            size_surface = 3
            size_fe = 4
            self.freedom = 3
        elif t == 6:
            self.fe_type = 'fe_2d_6'
            size_surface = 3
            size_fe = 6
            self.freedom = 2
        elif t == 8:
            self.fe_type = 'fe_3d_8'
            size_surface = 4
            size_fe = 8
            self.freedom = 3
        elif t == 10:
            self.fe_type = 'fe_3d_10'
            size_surface = 6
            size_fe = 10
            self.freedom = 3
        elif t == 24:
            self.fe_type = 'fe_2d_4'
            size_surface = 2
            size_fe = 4
            self.freedom = 2
        elif t == 34:
            self.fe_type = 'fe_1d_2'
            size_surface = 0
            size_fe = 2
            self.freedom = 1
        else:
            raise TFEMException('unknown_fe_err')
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









