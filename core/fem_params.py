#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#              Параметры решения задачи с помощью МКЭ
###################################################################

from core.fem_error import TException


# Тип задачи: статика или динамика
problem_type = [
    'static',
    'dynamic'
]

# Метод решения СЛАУ (точный, приближенный
solve_method = [
    'direct',
    'iterative'
]

# Стандартные имена функций (перемещения, деформации и напряжения) и их агрументов
std_name = [
    'x',    # 0  - идентификатор первого аргумента (x)
    'y',    # 1  - ... второго аргумента
    'z',    # 2  - ... третьего аргумента
    't',    # 3  - ... времени
    'U',    # 4  - компонента вектора перемещений по первому направлению (x)
    'V',    # 5  - ... по y
    'W',    # 6  - ... по z
    'Exx',  # 7  - компонента тензора нормальных деформаций по xx
    'Eyy',  # 8  - ... по yy
    'Ezz',  # 9  - ... по zz
    'Exy',  # 10 - компонента тензора тангенциальных деформаций по xу
    'Exz',  # 11 - ... по xz
    'Eyz',  # 12 - ... по yz
    'Sxx',  # 13 - компонента тензора нормальных напряжений по xx
    'Syy',  # 14 - ... по yy
    'Szz',  # 15 - ... по zz
    'Sxy',  # 16 - компонента тензора тангенциальных напряжений по xу
    'Sxz',  # 17 - ... по xz
    'Syz',  # 18 - ... по yz
    'Ut',   # 19 - скорость по x
    'Vt',   # 20 - ... по y
    'Wt',   # 21 - ... по z
    'Utt',  # 22 - ускорение по x
    'Vtt',  # 23 - ... по y
    'Wtt'   # 24 - ... по z
]

# Типы условий (краевых, упругих, ...)
type_condition = [
    'initial',
    'boundary',
    'volume',
    'surface',
    'concentrated',
    'pressure',
    'thickness',
    'young_modulus',
    'shear_modulus',
    'poisson_ratio',
    'temperature',
    'alpha',
    'density',
    'damping'
]


# Описание условия
class TParameter:
    def __init__(self):
        self.direct = 0         # Направление (номер функции, для которой задается условие): 0 - по x и т.д.
        self.expression = ''    # Функциональное выражение, определяющее значение условия (например, 10^5)
        self.predicate = ''     # Предикат отбора узлов
        self.type = ''          # Тип краевого условия


# Базовые параметры расчета задачи теории упругости с помощью МКЭ
class TFEMParams:
    def __init__(self):
        self.problem_type = ''  # Тип задачи
        self.solve_method = ''  # Метод решения СЛАУ
        self.width = 12         # Формат вывода результатов
        self.precision = 5
        self.eps = 1.0E-6       # Точность вычислений
        self.t0 = 0             # Начальный момент времени расчета
        self.t1 = 0             # Конечный момент времени расчета
        self.th = 0             # Шаг по времени
        self.names = std_name   # Список имен функций и их аргументов
        self.bc_list = []       # Список краевых условий
        self.var_list = {}      # Список вспомогательных переменных и их значений

    def __add_parameter(self, t, e, p, d=0):
        c = TParameter()
        c.type = t
        c.direct = d
        c.expression = e
        c.predicate = p
        self.bc_list.append(c)

    def add_boundary_condition(self, e, p, d):
        self.__add_parameter('boundary', e, p, d)

    def add_initial_condition(self, e, p, d):
        self.__add_parameter('initial', e, p, d)

    def add_volume_load(self, e, p, d):
        self.__add_parameter('volume', e, p, d)

    def add_surface_load(self, e, p, d):
        self.__add_parameter('surface', e, p, d)

    def add_concentrated_load(self, e, p, d):
        self.__add_parameter('concentrated', e, p, d)

    def add_pressure_load(self, e, p):
        self.__add_parameter('pressure', e, p, 0)

    def add_density(self, e, p):
        self.__add_parameter('density', e, p)

    def add_damping(self, e, p):
        self.__add_parameter('damping', e, p)

    def add_thickness(self, e, p):
        self.__add_parameter('thickness', e, p)

    def add_young_modulus(self, e, p):
        self.__add_parameter('young_modulus', e, p)

    def add_shear_modulus(self, e, p):
        self.__add_parameter('shear_modulus', e, p)

    def add_poisson_ratio(self, e, p):
        self.__add_parameter('poisson_ratio', e, p)

    def add_temperature(self, e, p):
        self.__add_parameter('temperature', e, p)

    def add_alpha(self, e, p):
        self.__add_parameter('alpha', e, p)

    def add_variable(self, var, val):
        self.var_list.setdefault(var, val)

    def check_params(self):
        if self.solve_method == '':
            raise TException('solve_method_err')
        if self.problem_type == '':
            raise TException('problem_type_err')
        is_young_modulus = False
        is_poisson_ratio = False
        for i in range(len(self.bc_list)):
            if self.bc_list[i].type == 'young_modulus':
                is_young_modulus = True
            if self.bc_list[i].type == 'poisson_ratio':
                is_poisson_ratio = True
        if not is_young_modulus or not is_poisson_ratio:
            raise TException('elasticity_err')
        if self.problem_type == 'dynamic' and (self.t0 == self.t1 or self.th <= 0):
            raise TException('time_err')
