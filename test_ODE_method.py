# Sam Alvis, January 2022

import time
from common_funcs import *
import copy


def main():
    # define body
    Earth = Body()
    Earth.mu = 3.9860E5
    Earth.R = 6371.0
    Earth.J2 = 1.9555E-3
    Earth.omega = 2 * math.pi / 88619.6
    Earth.max_atm_height = 200.0
    Earth.load_density_info("Atm_Data/0_ref.txt")

    # define craft
    SC = Craft()
    SC.beta_i = 50.0  # should be defined by ratio - half, etc. bigger the ratio, stronger the control event
    SC.beta_f = 150.0  # 140 by default?
    SC.t_eject = 0  # NOTE: here, this time will be adjusted and so is irrelevant

    # define x_0
    x_0 = [0, 0, 0, 0, 0, 0]

    method1 = 'RK45'  # Acceptable for single pass? Also could lower error tolerance in sim too?
    method2 = 'DOP853'

    num = 4
    for idx in range(num+1):
        # adjust x_0 slightly
        x_0 = [Earth.R + 200.0, 0, 0, 8, math.radians(-4 - 1*idx/num), math.pi/2]
        # run method 1
        tic1 = time.perf_counter()
        sol1 = orbit_sim(x_0, Earth, SC, run_time=1024, time_res=1, method=method1)
        toc1 = time.perf_counter()
        # run method 2
        tic2 = time.perf_counter()
        sol2 = orbit_sim(x_0, Earth, SC, run_time=1024, time_res=1, method=method2)
        toc2 = time.perf_counter()
        if sol1.hasImpacted:
            print("Method 1 Impacts!")
        if sol2.hasImpacted:
            print("Method 1 Impacts!")
        # compare times
        print(f"Method 1 ({method1}): {toc1-tic1} s; Method 2 ({method2}): {toc2-tic2} s")
        print(f"Method 1 ({method1}) is {toc2-tic2-toc1+tic1} s faster than Method 2 ({method2})")
       # compare results
        print(f"Method 1 ({method1}): {sol1.r[:,-1]} km; Method 2 ({method2}): {sol2.r[:,-1]} km")
        print(f"Error: {np.linalg.norm(sol1.r[:,-1] - sol2.r[:,-1])} km")


if __name__ == "__main__":
    main()
