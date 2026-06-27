import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.interpolate import make_interp_spline

strfuns = h5py.File('strfun.128.h5', 'r')
keys = list(strfuns.keys())
print(keys)
keys.remove('Dx')
Dx = np.array(strfuns['Dx'])
Dy = np.array(strfuns['Du_l3'])
x_data = np.log(Dx)
y_data = np.log(Dy)
# Create a spline interpolation of the data
spline = make_interp_spline(x_data, y_data, k=3)
first_derivative_spline = spline.derivative(nu=1)
x = Dx
exp = first_derivative_spline(np.log(Dx))
plt.plot(x, exp, label='Num', color='blue')
plt.plot(x, np.ones_like(x), label='exp=1', color='red', linestyle='--')
plt.xlabel('dlog(Dx)')
plt.ylabel('dlog(Du_l3)')
plt.title('Spline Interpolation of log(Du_l3) vs log(Dx)')
plt.legend()
plt.grid()
plt.show()