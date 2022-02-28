# Sam Alvis, January 2022

import time
from common_funcs import *
import copy


def main():
    # define body
    Earth_0 = Body()
    Earth_0.mu = 3.9860E5
    Earth_0.R = 6371.0
    Earth_0.J2 = 1.9555E-3
    Earth_0.omega = 2 * math.pi / 88619.6
    Earth_0.max_atm_height = 200.0
    Earth_0.load_density_info("Atm_Data/0_ref.txt")

    # define craft
    SC_0 = Craft()
    SC_0.beta_i = 50.0  # should be defined by ratio - half, etc. bigger the ratio, stronger the control event
    SC_0.beta_f = 150.0  # 140 by default?
    SC_0.t_eject = 0  # NOTE: here, this time will be adjusted and so is irrelevant

    # define params
    a_f_targ = Earth_0.R + 1000.0  # km

    # define x_0
    x_0 = [Earth_0.R + 200.0, 0, 0, 8, math.radians(-4.5), math.pi/2]

    # define target altitude and run function, then plot results
    Earth_0.scale_factor = 1
    tic = time.perf_counter()
    solution = eject_calc(a_f_targ, x_0, Earth_0, SC_0)
    toc = time.perf_counter()
    print(f"t_eject = {solution.craft_obj.t_eject} [took {toc-tic:0.4f} seconds]")
    plot_case(solution)


def main2():
    # define body
    Earth_0 = Body()
    Earth_0.mu = 3.9860E5
    Earth_0.R = 6371.0
    Earth_0.J2 = 1.9555E-3
    Earth_0.omega = 2 * math.pi / 88619.6
    Earth_0.max_atm_height = 200.0
    Earth_0.load_density_info("Atm_Data/0_ref.txt")
    Earth_0.scale_factor = 1

    # define craft
    SC_0 = Craft()
    SC_0.beta_i = 50.0  # should be defined by ratio - half, etc. bigger the ratio, stronger the control event
    SC_0.beta_f = 150.0  # 140 by default?
    SC_0.t_eject = 0  # NOTE: here, this time will be adjusted and so is irrelevant

    # define params
    a_f_targ = Earth_0.R + 1000.0  # km

    # define x_0
    x_0 = [Earth_0.R + 200.0, 0, 0, 8, math.radians(-4.5), math.pi/2]

    t_arr = np.linspace(0, 1024, 5)
    solutions = []
    tic = time.perf_counter()
    for t_ej in t_arr:
        print(t_ej)
        SC_0.t_eject = t_ej
        tic2 = time.perf_counter()
        this_sol = orbit_sim(x_0, Earth_0, copy.copy(SC_0), run_time=1024, time_res=1)
        toc2 = time.perf_counter()
        tic3 = time.perf_counter()
        solutions = np.append(solutions, this_sol)
        toc3 = time.perf_counter()
        print(f"{toc2-tic2}; {toc3-tic3}; {this_sol.craft_obj.t_eject}")

    toc = time.perf_counter()
    print(f"[took {toc-tic:0.4f} seconds]")

    t_ej_arr = [sol.craft_obj.t_eject for sol in solutions]
    a_f_arr = [sol.a_t[-1]-a_f_targ for sol in solutions]

    plt.plot(t_ej_arr, a_f_arr)
    plt.show()


if __name__ == "__main__":
    main()
