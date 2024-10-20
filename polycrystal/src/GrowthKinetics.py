import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

data = pd.read_csv("./AvgGrainSize.csv")

grain_size = data.iloc[:,1]
t = data.iloc[:,0]

R_0 = grain_size[0]

def growth_rate(x_range, param1, param2):
    
    y = (R_0**param1 + param2 * x_range)**(1/param1)

    return y

#popt, pcov = sp.optimize.curve_fit(growth_rate, t , grain_size, [3 , 150])

plt.scatter(t , grain_size, color='green', s = 0.5)
#plt.scatter(t, growth_rate(t, 3, 150), color='red', s = 0.5)
plt.show()

