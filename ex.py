from chain import ConstName, Var
from formatter_lsode import LsodeFormatter

y1 = Var('y1')
y2 = Var('y2')
y3 = Var('y3')

c1 = ConstName('c1')
c2 = ConstName('c2')
c3 = ConstName('c3')

lf = LsodeFormatter((y1, y2, y3), (
    -c1*y1 + c2*y2*y3,
    c1*y1 - c2*y2*y3 - c3*y2**2,
    c3*y2**2,
    ))

with open('fex.f90', 'w') as f:
    f.write(lf.make_functions_file())
with open('ex.f90', 'w') as f:
    f.write(lf.make_program_file(
                y0=[1.0, .0, .0],
                consts={
                    c1: .04,
                    c2: 1e4,
                    c3: 3e7,
                    },
                tout=.4,
                rtol=1e-4,
                atol=[1e-6, 1e-10, 1e-6],
                num_steps=12,
                ))
