#%% 

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# %%

x = np.linspace(0,100)
y = [xi + np.random.normal(0, 0.1) for xi in x]

fig = plt.figure()
ax = fig.gca()

ax.plot(x, x, color = "black")
ax.scatter(x, y, color = "black")


sns.despine()



plt.plot(y)