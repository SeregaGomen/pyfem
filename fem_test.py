#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fem_defs import DIR_X, DIR_Y, DIR_Z
from fem_error import TFEMException
from fem_object import TObject


def body1d():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/body1d.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
    obj.add_boundary_condition('0', 'x=0', DIR_X)
    obj.add_volume_load('-1.0E+5', '', DIR_X)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()

def cube():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/cube.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
    obj.add_boundary_condition('0', 'z=0', DIR_X | DIR_Y | DIR_Z)
    obj.add_volume_load('-1.0E+5', '', DIR_Z)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()


def beam():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/beam.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
    obj.add_boundary_condition('0', 'y=0', DIR_X | DIR_Y | DIR_Z)
    obj.add_volume_load('-1.0E+5', '', DIR_Y)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()


def console():
    obj = TObject()
    e = [6.5E+10]
    m = [0.3]
    obj.set_mesh('mesh/console.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
    obj.add_boundary_condition('0', 'x=0', DIR_X | DIR_Y)
    obj.add_volume_load('-1.0E+5', '', DIR_X)
    if obj.calc():
        obj.calc_results()
        obj.set_width(10)
        obj.set_precision(5)
        obj.print_result()

try:
    # body1d()
    # cube()
    console()
except TFEMException as err:
    err.print_error()
