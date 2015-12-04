#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fem_defs import DIR_X, DIR_Y, DIR_Z
from fem_parser import TParser
from fem_error import TFEMException
from fem_object import TObject

# code = '((x^2+ + y^2 <= R^2) or (x >= 0) or (x <= y))'

# p = TParser()


# p.add_variable('x', 0)
# p.add_variable('y', 0)
# p.add_variable('R', 2)
# p.set_code(code)
# if p.error == '':
#    print(p.run())


obj = TObject()
e = [6.5E+10]
m = [0.3]
try:
#    obj.set_mesh('body1d.trpa')
    obj.set_mesh('body.trpa')
    obj.set_problem_type('static')
    obj.set_solve_method('direct')
    obj.set_elasticity(e, m)
#    obj.add_boundary_condition('0', 'x=0', DIR_X)
#    obj.add_volume_condition('1.0E+6', '', DIR_X)
    obj.add_boundary_condition('0', 'y=0', DIR_X | DIR_Y)
    obj.add_volume_condition('1.0E+6', '', DIR_Y)
    obj.calc()

#    for i in range(0,len(obj.__global_matrix__[0])):
#        print(obj.__global_matrix__[i])
#    print(obj.__global_vector__)


except TFEMException as e:
    e.print_error()



