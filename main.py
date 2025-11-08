import math
from chain import (
    Arctan, ConstName, Cos, Sin, Sq, Sqrt, NamedFunction, ExplicitFunction,
    Var, PiecewiseFunction, ZERO
)

from formatter_lsode import LsodeFormatter

PI = ConstName('PI')

XG = Var('XG')
YG = Var('YG')
a = Var('a')

XG_dot = Var('XG_dot')
YG_dot = Var('YG_dot')
a_dot = Var('a_dot')

xG = ConstName('xG')
xT = ConstName('xT')

XTL = ExplicitFunction('XTL', XG + (xG - xT)*Cos(a))
YTL = ExplicitFunction('YTL', YG - (xG - xT)*Sin(a))

LTL = ExplicitFunction('LTL', Sqrt(Sq(XTL) + Sq(YTL)))
kTL = ConstName('kTL')
LTL0 = ConstName('LTL0')
FTL = ExplicitFunction('FTL', (LTL - LTL0)*kTL)

tang = ExplicitFunction('tang', YTL/XTL)
cosg = ExplicitFunction('cosg', XTL/LTL)
sing = ExplicitFunction('sing', YTL/LTL)
gam = ExplicitFunction('gamma', Arctan(tang))

FPDmax = ConstName('FPDmax')
FPDbow = ConstName('FPDbow')
FPDaft = ConstName('FPDaft')
FPD1 = ExplicitFunction('FPD1', FPDmax*(FPDbow + (1-FPDbow)*Sq(Sin(a))))
FPD2 = ExplicitFunction('FPD2', FPDmax*(FPDaft + (1-FPDaft)*Sq(Sin(a))))

FPD = ExplicitFunction('FPD',
        PiecewiseFunction.factory({
            a << PI/2: FPD1,
            a << PI:   FPD2,
        }, ZERO))

xCPmax = ConstName('xCPmax')
xCPmed = ConstName('xCPmed')
xCP = ExplicitFunction('xCP', xCPmax*Sq(Cos(a/2)) + (xCPmed - xCPmax/2)*Sq(Sin(a)))

FPLmax = ConstName('FPLmax')
FPL1 = ExplicitFunction('FPL1', FPLmax*Sq(Sin(2*a)))
FPL2 = ExplicitFunction('FPL2', FPLmax*Sq(Cos(4*(a-PI/4))))
FPL = ExplicitFunction('FPL',
       PiecewiseFunction.factory({
        a << PI/4:   FPL1,
        a << 3*PI/8: FPL2,
        }, ZERO))

RX = ExplicitFunction('RX', FPD - FTL*cosg)
RY = ExplicitFunction('RY', FPL + FTL*sing)
M = ExplicitFunction('M', FTL*Sin(a + gam)*(xG - xT)   \
                        - FPD*(xG - xCP)*Sin(a)     \
                        + FPL*(xG - xCP)*Cos(a))

m = ConstName('m')
Iz = ConstName('Iz')

ode_y = [XG, YG, a, XG_dot, YG_dot, a_dot]
ode_f = [XG_dot, YG_dot, a_dot, RX/m, RY/m, M/Iz]

LWL = 10.
tf_to_N = 9.80665e3

# Determination of kTL
#
# Not easy, but I found a graph in TUP4 that will help: 
#   Fig. 7.27: load (% of breaking load) vs deformation (%) - TUP4 (p. 193)
# Let 'f' be the ratio (derivative) of that graph.
#   k = EA/L, but will be approximated by f*Fmax/L, so EA = f*Fmax
# I will determine Fmax from parts of TUP4 (for some given reasonable diameter around 40mm)
# Below we find EAmin (nylon) and EAmax (steel), and use a % in our simulation.
EAmin = (.65/.33) * 250e3
EAmax = (.65/.02) * 110. * tf_to_N
# EAmax/EAmin ~= 71
# Nylon is too stretchy, so we'll be like 60% steel

_LTL0 = 100.

consts = {PI: math.pi,
          xG: .5 * LWL,
          xT: .3 * LWL,
          kTL: (EAmin + .60*(EAmax-EAmin))/_LTL0,
          LTL0: _LTL0,
          xCPmax: .5*LWL,
          xCPmed: .3*LWL,
          FPLmax: 5.*tf_to_N,
          FPDmax: 10.*tf_to_N,
          FPDbow: .1,
          FPDaft: .2,
          m: 1e4,
          Iz: 1e5,
          FPDmax: 10. * tf_to_N,
          FPLmax: 3. * tf_to_N,
          }

lf = LsodeFormatter(ode_y, ode_f)

with open('fmain.f90', 'w') as f:
    f.write(lf.make_functions_file())
with open('main.f90', 'w') as f:
    f.write(lf.make_program_file(
                y0=[100.0, .0, .0, .0, .0, .0],
                consts=consts,
                tout=600.,
                rtol=1e-4,
                atol=[1., 1., .01, 1., 1., .01],
                ))
