from chain import Const, Var
from formatter_lsode import LsodeFormatter

y1 = Var('y1')
y2 = Var('y2')
y3 = Var('y3')

c1 = Const(.04)
c2 = Const(1e4)
c3 = Const(3e7)

lf = LsodeFormatter((y1, y2, y3), (
    -c1*y1 + c2*y2*y3,
    c1*y1 - c2*y2*y3 - c3*y2**2,
    c3*y2**2,
    ))

with open('fex.f90', 'w') as f:
    f.write(lf.make_functions_file())
with open('ex.f90', 'w') as f:
    f.write(lf.make_main(['1.0_dp','.0_dp','.0_dp'],
                       '.0_dp',
                       '.4_dp',
                       '1e-4_dp',
                       ['1e-6_dp','1e-10_dp','1e-6_dp']))
