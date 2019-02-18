#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#           Конечно-элементная модель объекта расчета
###################################################################

from core.fem_error import TFEMException


# Типы конечных элементов
FEType = [
    'fe_1d_2',
    'fe_2d_3',
    'fe_2d_6',
    'fe_2d_3_p',
    'fe_2d_6_p',
    'fe_2d_3_s',
    'fe_2d_6_s',
    'fe_2d_4',
    'fe_2d_4_p',
    'fe_3d_4',
    'fe_3d_8',
    'fe_3d_10'
]


class TMesh:
    def __init__(self):
        self.mesh_file = ''         # Имя файла с данными о геометрической модели
        self.fe_type = ''           # Тип КЭ
        self.x = []                 # Координаты узлов
        self.be = []                # Связи граничных элементов
        self.fe = []                # Связи в КЭ
        self.freedom = 0            # Кол-во степеней свободы
        self.dimension = 0          # Размерность КЭ
        self.surface_links = []     # Список граничных элементов, соответствующих каждому узлу

    @staticmethod
    def get_fe_data(t):
        if t == 3:
            return 'fe_2d_3', 2, 3, 2, 2
        elif t == 4:
            return 'fe_3d_4', 3, 4, 3, 3
        elif t == 6:
            return 'fe_2d_6', 3, 6, 2, 2
        elif t == 8:
            return 'fe_3d_8', 4, 8, 3, 3
        elif t == 10:
            return 'fe_3d_10', 6, 10, 3, 3
        elif t == 24:
            return 'fe_2d_4', 2, 4, 2, 2
        elif t == 34:
            return 'fe_1d_2', 1, 2, 1, 1
        elif t == 123:
            return 'fe_2d_3_p', 0, 3, 3, 2
        elif t == 124:
            return 'fe_2d_4_p', 0, 4, 3, 2
        elif t == 125:
            return 'fe_2d_6_p', 0, 6, 3, 2
        elif t == 223:
            return 'fe_2d_3_s', 0, 3, 6, 3
        elif t == 224:
            return 'fe_2d_4_s', 0, 4, 6, 3
        elif t == 225:
            return 'fe_2d_6_s', 0, 6, 6, 3
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
        self.surface_links = [[] for i in range(n)]
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
                k = int(lines[index].split()[j])
                row.append(k)
                self.surface_links[k].append(i)
            self.be.append(row)
            index += 1
        if self.fe_type == 'fe_2d_3_p' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_2d_3_s' or \
                self.fe_type == 'fe_2d_4_p' or self.fe_type == 'fe_2d_6_s' or self.fe_type == 'fe_2d_4_s':
            self.be = self.fe

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
        elif self.fe_type == 'fe_2d_3_s':
            return 'triangular shell element (3 nodes)'
        elif self.fe_type == 'fe_2d_6_s':
            return 'triangular shell element (6 nodes)'
        elif self.fe_type == 'fe_2d_4_p':
            return 'plate element (4 nodes)'
        elif self.fe_type == 'fe_2d_4_s':
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
        return True if self.fe_type == 'fe_2d_3_s' or self.fe_type == 'fe_2d_4_s' or self.fe_type == 'fe_2d_6_s' \
            else False

    def base_be_size(self):
        if self.fe_type == 'fe_2d_6':
            return 2
        elif self.fe_type == 'fe_3d_10' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_2d_6_s':
            return 3
        return len(self.fe[0])

    def base_fe_size(self):
        if self.fe_type == 'fe_2d_6' or self.fe_type == 'fe_2d_6_p' or self.fe_type == 'fe_2d_6_s':
            return 3
        elif self.fe_type == 'fe_3d_10':
            return 4
        return len(self.fe[0])

    # Процедура поиска граничных ребер (граней) у КЭ
    def get_surface_list(self, index):
        s_list = []
        if self.fe_type == 'fe_1d_2':
            s_index = [[0], [1]]
        elif self.fe_type == 'fe_2d_3' or self.fe_type == 'fe_2d_6':
            s_index = [[0, 1], [1, 2], [2, 0]]
        elif self.fe_type == 'fe_2d_4':
            s_index = [[0, 1], [1, 2], [2, 3], [3, 0]]
        elif self.fe_type == 'fe_3d_4' or self.fe_type == 'fe_3d_10':
            s_index = [[0, 1, 3], [1, 2, 3], [2, 0, 3], [0, 1, 2]]
        elif self.fe_type == 'fe_3d_8':
            s_index = [[0, 1, 2, 3], [4, 5, 6, 7], [1, 2, 6, 5], [3, 0, 4, 7], [0, 1, 5, 4], [2, 3, 7, 6]]
        else:
            s_index = []
        # Формируем список граничных элементов, вершинами которых являются узлы треугольника
        be_list = []
        for i in range(len(self.fe[index])):
            for j in range(len(self.surface_links[self.fe[index][i]])):
                be_list.append(set(self.be[self.surface_links[self.fe[index][i]][j]]))
        # Проверяем, есть ли ребра (грани) КЭ в этом списке
        for i in range(len(s_index)):
            surface = set([])
            for j in range(len(s_index[i])):
                surface.add(self.fe[index][s_index[i][j]])
            for j in range(len(be_list)):
                if surface == be_list[j]:
                    s_list.append(i)
                    break
        if (len(s_list)):
            print(s_list)
        return s_list
