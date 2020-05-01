# visualize_google_matrix.py
# honestly nothing about this is specific to Google matrices it just makes a
# heatmap

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

data = np.load('sars-covid_google_mat.npy')
sns.heatmap(data)
plt.show()
