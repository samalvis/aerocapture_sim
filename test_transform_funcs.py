# Sam Alvis, October 2021
import random

from Archived.general_funcs import *


# define body
Test = Body()
Test.mu = 4.2828E4
Test.R = 3396.2
Test.J2 = 1.9555E-3
Test.omega = 2 * math.pi / 88619.6
Test.max_atm_height = 100
Test.load_density_info("mars_atm_data.csv")

# define craft
SC = Craft()
SC.beta = 140

x_test_inert = [3400, 0, 0, 0, 0, 1]
x_test = car_2_loc(x_test_inert, Test, 100)
x_test_fin = loc_2_car(x_test, Test, 100)
print(x_test_inert)
print(x_test)
print(x_test_fin)

max_num = 100
X_i = np.zeros((6, max_num))
X_m = np.zeros((6, max_num))
X_f = np.zeros((6, max_num))
err_x = np.zeros(max_num)
err_v = np.zeros(max_num)

t_test = 0

for num in range(max_num):
    X_i[:,num] = [random.uniform(-10000,10000), random.uniform(-10000,10000), random.uniform(-10000,10000), random.uniform(-10,10), random.uniform(-10,10), random.uniform(-10,10)]
    X_m[:,num] = car_2_loc(X_i[:,num], Test, t_test)
    X_f[:,num] = loc_2_car(X_m[:,num], Test, t_test)
    err_x[num] = np.linalg.norm(X_i[0:3,num] - X_f[0:3,num])
    err_v[num] = np.linalg.norm(X_i[3:6,num] - X_f[3:6,num])

#print(np.transpose(X_i[:,0:1]))
#print(np.transpose(X_m[:,0:1]))
#print(np.transpose(X_f[:,0:1]))
print(max(err_x))
print(max(err_v))


