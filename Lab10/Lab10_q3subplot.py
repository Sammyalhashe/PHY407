#!/usr/bin/env python
import matplotlib.pyplot as plt

from Lab10_q3b import Iarray1
from Lab10_q3c import Iarray2

graph, (ax1, ax2) = plt.subplots(1,2)

ax1.hist(Iarray1,10,range=[0.8, 0.88])
ax1.set_xlabel('Mean Value Method Result')

ax2.hist(Iarray2,10,range=[0.8, 0.88])
ax2.set_xlabel('Importance Sampling Result')

graph.text(0.32, 0.98, '100 Iterations with Different Methods', va='center')

plt.tight_layout()
plt.show()