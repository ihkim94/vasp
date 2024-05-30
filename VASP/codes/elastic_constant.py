import numpy as np
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline

x = [0.001, -0.001, 0, 4, 5, 8, 10]
y = [-41.629630, 3.8, -41.629630, 16, 24, 65, 99.2]

model = make_pipeline(PolynomialFeatures(2), LinearRegression())
model.fit(np.array(x).reshape(-1, 1), y)
x_reg = np.arange(11)
y_reg = model.predict(x_reg.reshape(-1, 1))

plt.scatter(x, y)
plt.plot(x_reg, y_reg)
plt.show()