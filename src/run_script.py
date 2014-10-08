import geometry
import MOC
g = MOC.makeGeometry()
m = MOC.MOC(g)
m.makeTracks()
m.makeReflective()
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
for i in m.tracks3D:
    for p in i:
        for a in i:
            for d in a:
                for t in d:
                    ax.plot([t.r_in.x, t.r_out.x], [t.r_in.y, t.r_out.y], [t.r_in.z, t.r_out.z], 'k')

plt.show()

