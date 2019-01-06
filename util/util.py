#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#                    Вспомогательные процедуры
#######################################################################

import math


# Создание сетки для квадратной пластины
def create_plate_mesh_4():
    x_min = [-0.5, -0.5]
    x_max = [0.5, 0.5]
    n = 200
    h = [(x_max[0] - x_min[0])/n, (x_max[1] - x_min[1])/n]
    index = []
    x = []
    counter = 0
    for i in range(0, n + 1):
        c_index = []
        for j in range(0, n + 1):
            c_index.append(counter)
            counter += 1
            x.append([x_min[0] + i*h[0], x_min[1] + j*h[1]])
        index.append(c_index)
    with open('mesh/plate4-1.0.trpa', 'w') as file:
        file.write('124\n')
        file.write(str(counter) + '\n')
        for i in range(0, len(x)):
            file.write(str(x[i][0]) + ' ' + str(x[i][1]) + '\n')
        file.write(str(n**2) + '\n')
        for i in range(0, len(index) - 1):
            for j in range(0, len(index) - 1):
                file.write(str(index[i][j]) + ' ' + str(index[i][j + 1]) + ' ' + str(index[i + 1][j + 1]) + ' ' +
                           str(index[i + 1][j]) + '\n')
        file.write('0\n')
    return


# Создание сетки для оболочки в форме полого цилиндра
def create_shell_mesh_4():
    r = 3.98/2
    height = 4.014
    n_xy = 100
    n_z = 50
    d_fi = 2*math.pi/n_xy
    d_h = height/n_z
    index = []
    x = []
    counter = 0
    for i in range(0, n_z + 1):
        c_index = []
        for j in range(0, n_xy):
            c_index.append(counter)
            counter += 1
            x.append([r*math.cos(j*d_fi), r*math.sin(j*d_fi), i*d_h])
        index.append(c_index)
    with open('mesh/shell4-1.0.trpa', 'w') as file:
        file.write('224\n')
        file.write(str(counter) + '\n')
        for i in range(0, len(x)):
            file.write(str(x[i][0]) + ' ' + str(x[i][1]) + ' ' + str(x[i][2]) + '\n')
        file.write(str(n_xy*n_z) + '\n')
        for i in range(0, len(index) - 1):
            for j in range(0, len(index[0]) - 1):
                file.write(str(index[i][j]) + ' ' + str(index[i][j + 1]) + ' ' + str(index[i + 1][j + 1]) + ' ' +
                           str(index[i + 1][j]) + '\n')
            file.write(str(index[i][j + 1]) + ' ' + str(index[i][0]) + ' ' + str(index[i + 1][0]) + ' ' +
                       str(index[i + 1][j + 1]) + '\n')
        file.write('0\n')
    return


# Конвертация файла данных gmsh в формат trpa
def convert_msh_2_trpa(file_msh, file_trpa):
    with open(file_msh, 'r') as file:
        line = file.readline()
        if line != '$MeshFormat\n':
            print('Error file format %s' % file_msh)
            return False
        while line != '':
            line = file.readline()
            if line == '$EndMeshFormat\n':
                break
        line = file.readline()
        if line == '$Entities\n':
            while line != '':
                line = file.readline()
                if line == '$EndEntities\n':
                    break
        x = []
        t_index = []
        line = file.readline()
        if line == '$Nodes\n':
            line = file.readline()
            while line != '':
                line = file.readline()
                if line == '$EndNodes\n':
                    break
                sx = line.split(' ')
                n = int(sx[3])  # Количество координат в текущем блоке
                for i in range(0, n):
                    line = file.readline()
                    sx = line.split(' ')
                    cx = []
                    t_index.append(int(sx[0]))
                    for j in range(1, len(sx)):
                        cx.append(float(sx[j]))
                    x.append(cx)
        # Переиндексируем номера узлов
        max_index = t_index[len(t_index) - 1]
        index = [0]*max_index
        for i in range(len(t_index)):
            index[t_index[i] - 1] = i
        fe = []
        be = []
        line = file.readline()
        if line == '$Elements\n':
            line = file.readline()
            while line != '':
                line = file.readline()
                if line == '$EndElements\n':
                    break
                sl = line.split(' ')
                n = int(sl[3])  # Количество элементов в текущем блоке
                for i in range(0, n):
                    line = file.readline().replace('\n', '').strip()
                    se = line.split(' ')
                    if len(se) == 2:
                        continue
                    if len(se) == 3:    # ГЭ в форме отрезка
                        cbe = []
                        for j in range(1, len(se)):
                            cbe.append(index[int(se[j]) - 1])
                        be.append(cbe)
                    if len(se) == 4 or len(se) == 5:    # КЭ в форме треугольника или четрыехугольника
                        cfe = []
                        for j in range(1, len(se)):
                            cfe.append(index[int(se[j]) - 1])
                        fe.append(cfe)
    with open(file_trpa, 'w') as file:
        file.write('24\n')
        file.write('%d\n' % len(x))
        for i in range(0, len(x)):
            for j in range(0, len(x[i]) - 1):
                file.write(str(x[i][j]) + ' ')
            file.write('\n')
        file.write('%d\n' % len(fe))
        for i in range(0, len(fe)):
            for j in range(0, len(fe[i])):
                file.write(str(fe[i][j]) + ' ')
            file.write('\n')
        file.write('%d\n' % len(be))
        for i in range(0, len(be)):
            for j in range(0, len(be[i])):
                file.write(str(be[i][j]) + ' ')
            file.write('\n')
    return True


# Вспомогательная процедура поиска ребра в списке обработанных
def is_find(le, e):
    for i in range(0, len(le)):
        if le[i][0] != e[0]:
            continue
        if le[i][1] == e[1]:
            return True, le[i][2]
    return False, -1


# Конвертация треугольных КЭ в квадратичные
def mesh_convert_2d_3_2_6(file_linear, file_quadric):
    x = []
    be = []
    fe = []
    # Считываем исходную сетку
    with open(file_linear, 'r') as file:
        line = file.readline()
        if int(line) != 3:
            print('Mesh error: incorrect FE type')
            return False
        # Считываем узлы
        n = int(file.readline())
        for i in range(0, n):
            line = file.readline().replace('\n', '').strip()
            s = line.split(' ')
            x.append([float(s[0]), float(s[1])])
        # Считываем КЭ
        n = int(file.readline())
        for i in range(0, n):
            line = file.readline().replace('\n', '').strip()
            s = line.split(' ')
            fe.append([int(s[0]), int(s[1]), int(s[2])])
        # Считываем ГЭ
        n = int(file.readline())
        for i in range(0, n):
            line = file.readline().replace('\n', '').strip()
            s = line.split(' ')
            be.append([int(s[0]), int(s[1])])
    # Конвертируем...
    n_fe = []
    n_be = []
    edge = []
    # Преобразуем границу
    num = len(x)
    for i in range(0, len(be)):
        xp = [(x[be[i][0]][0] + x[be[i][1]][0]) / 2, (x[be[i][0]][1] + x[be[i][1]][1]) / 2]
        x.append(xp)
        n_be.append([be[i][0], be[i][1], num])
        e = sorted(be[i])
        edge.append([e[0], e[1], num])
        num += 1
    # Преобразуем КЭ
    for i in range(0, len(fe)):
        index = [[fe[i][0], fe[i][1]], [fe[i][1], fe[i][2]], [fe[i][2], fe[i][0]]]
        add_fe = []
        for j in range(0, 3):
            ret, c_index = is_find(edge, sorted(index[j]))
            if ret is not True:
                xp = [(x[index[j][0]][0] + x[index[j][1]][0]) / 2, (x[index[j][0]][1] + x[index[j][1]][1]) / 2]
                x.append(xp)
                add_fe.append(num)
                e = sorted(index[j])
                edge.append([e[0], e[1], num])
                num += 1
            else:
                add_fe.append(c_index)
        n_fe.append([fe[i][0], fe[i][1], fe[i][2], add_fe[0], add_fe[1], add_fe[2]])
    # Выводим результат
    with open(file_quadric, 'w') as file:
        file.write('6\n')
        file.write('%d\n' % len(x))
        for i in range(0, len(x)):
            for j in range(0, len(x[i])):
                file.write(str(x[i][j]) + ' ')
            file.write('\n')
        file.write('%d\n' % len(n_fe))
        for i in range(0, len(n_fe)):
            for j in range(0, len(n_fe[i])):
                file.write(str(n_fe[i][j]) + ' ')
            file.write('\n')
        file.write('%d\n' % len(n_be))
        for i in range(0, len(n_be)):
            for j in range(0, len(n_be[i])):
                file.write(str(n_be[i][j]) + ' ')
            file.write('\n')
    return True


# Конвертация линейного тетраэдра в квадратичный
def mesh_convert_3d_4_2_10(file_linear, file_quadric):
    x = []
    be = []
    fe = []
    # Считываем исходную сетку
    with open(file_linear, 'r') as file:
        line = file.readline()
        if int(line) != 4:
            print('Mesh error: incorrect FE type')
            return False
        # Считываем узлы
        n = int(file.readline())
        for i in range(0, n):
            line = file.readline().replace('\n', '').strip()
            s = line.split(' ')
            x.append([float(s[0]), float(s[1]), float(s[2])])
        # Считываем КЭ
        n = int(file.readline())
        for i in range(0, n):
            line = file.readline().replace('\n', '').strip()
            s = line.split(' ')
            fe.append([int(s[0]), int(s[1]), int(s[2]), int(s[3])])
        # Считываем ГЭ
        n = int(file.readline())
        for i in range(0, n):
            line = file.readline().replace('\n', '').strip()
            s = line.split(' ')
            be.append([int(s[0]), int(s[1]), int(s[2])])
    # Конвертируем...
    n_fe = []
    n_be = []
    edge = []
    # Преобразуем границу
    num = len(x)
    for i in range(0, len(be)):
        index = [[be[i][0], be[i][1]], [be[i][1], be[i][2]], [be[i][0], be[i][2]]]
        a_be = []
        for j in range(0, 3):
            ret, c_index = is_find(edge, sorted(index[j]))
            if ret is not True:
                xp = [(x[index[j][0]][0] + x[index[j][1]][0]) / 2, (x[index[j][0]][1] + x[index[j][1]][1]) / 2,
                      (x[index[j][0]][2] + x[index[j][1]][2]) / 2]
                x.append(xp)
                a_be.append(num)
                e = sorted(index[j])
                edge.append([e[0], e[1], num])
                num += 1
            else:
                a_be.append(c_index)
        n_be.append([be[i][0], be[i][1], be[i][2], a_be[0], a_be[1], a_be[2]])

    # Преобразуем КЭ
    for i in range(0, len(fe)):
        index = [[fe[i][0], fe[i][1]], [fe[i][1], fe[i][2]], [fe[i][0], fe[i][2]], [fe[i][2], fe[i][3]],
                 [fe[i][1], fe[i][3]], [fe[i][0], fe[i][3]]]
        a_fe = []
        for j in range(0, 6):
            ret, c_index = is_find(edge, sorted(index[j]))
            if ret is not True:
                xp = [(x[index[j][0]][0] + x[index[j][1]][0]) / 2, (x[index[j][0]][1] + x[index[j][1]][1]) / 2,
                      (x[index[j][0]][2] + x[index[j][1]][2]) / 2]
                x.append(xp)
                a_fe.append(num)
                e = sorted(index[j])
                edge.append([e[0], e[1], num])
                num += 1
            else:
                a_fe.append(c_index)
        n_fe.append([fe[i][0], fe[i][1], fe[i][2], fe[i][3], a_fe[0], a_fe[1], a_fe[2], a_fe[3], a_fe[4], a_fe[5]])
    # Выводим результат
    with open(file_quadric, 'w') as file:
        file.write('10\n')
        file.write('%d\n' % len(x))
        for i in range(0, len(x)):
            for j in range(0, len(x[i])):
                file.write(str(x[i][j]) + ' ')
            file.write('\n')
        file.write('%d\n' % len(n_fe))
        for i in range(0, len(n_fe)):
            for j in range(0, len(n_fe[i])):
                file.write(str(n_fe[i][j]) + ' ')
            file.write('\n')
        file.write('%d\n' % len(n_be))
        for i in range(0, len(n_be)):
            for j in range(0, len(n_be[i])):
                file.write(str(n_be[i][j]) + ' ')
            file.write('\n')
    return True


# convert_msh_2_trpa('/home/serg/work/Qt/QFEM/QFEM/mesh/tank-new/gmsh/quad-1.msh', 'mesh/quad-4-1.trpa')
convert_msh_2_trpa('../mesh/quad-1.msh', '../mesh/quad-4.trpa')
# create_shell_mesh_4()
# create_plate_mesh_4()
# mesh_convert_2d_3_2_6('../mesh/quad-3.trpa', '../mesh/quad-6.trpa')

# mesh_convert_2d_3_2_6('../mesh/console.trpa', '../mesh/console-6.trpa')
# mesh_convert_3d_4_2_10('../mesh/beam.trpa', '../mesh/beam-10.trpa')
