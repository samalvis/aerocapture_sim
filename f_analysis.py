# Sam Alvis, January 2022

from common_funcs import *
import copy


def main():
    """Description"""

    # define body
    Earth_0 = Body()
    Earth_0.mu = 3.9860E5
    Earth_0.R = 6371.0
    Earth_0.J2 = 1.9555E-3
    Earth_0.omega = 2 * math.pi / 88619.6
    Earth_0.max_atm_height = 200.0
    Earth_0.load_density_info("Atm_Data/0_ref.txt")  # TODO: fix density function itself! Not log!

    # define craft
    SC_0 = Craft()
    SC_0.beta_i = 50.0  # should be defined by ratio - half, etc. bigger the ratio, stronger the control event
    SC_0.beta_f = 150.0  # 140 by default?
    SC_0.t_eject = 0  # NOTE: here, this time will be adjusted and so is irrelevant

    # define params
    a_f_targ = Earth_0.R + 1000.0  # km

    # define x_0
    x_0 = [Earth_0.R + 200.0, 0, 0, 8, math.radians(-4.5), math.pi/2]

    print("Starting f_list generation")

    optimal_solution_arr = []
    f = 1
    f_inc = 0.05
    # increase f up until crash  # TODO: is this count up count down method worth it? It's slighlty more efficient but could have errors if it starts out of the capture range...
    while f <= 1.1:  # NOTE: (removed and still works?) this 2 is arbitrary, I just want it to be bounded
        print(f)
        # adjust body
        Earth = copy.copy(Earth_0)
        SC = copy.copy(SC_0)
        # find t_ej, a_star(t), a_dot_star(t) (?)
        Earth.scale_factor = f
        solution = eject_calc_2(a_f_targ, x_0, Earth, SC)
        optimal_solution_arr.append(solution)  # TODO: fix hard code here!
        print(solution.craft_obj.t_eject)
        if solution.craft_obj.t_eject == 0:  # NOTE: i don't know if this is working - should it be eject time? final apoapsis?
            break                 # b/c it's still capturing, just not reaching the high enough orbit
        f = f + f_inc

    print("Halfway done with f_list generation")

    f = 1 - f_inc
    # decrease f up until escape or 0 atm
    while f >= 0.9:
        print(f)
        # adjust body
        Earth = copy.copy(Earth_0)
        SC = copy.copy(SC_0)
        # find t_ej, a_star(t), a_dot_star(t) (?)
        Earth.scale_factor = f
        solution = eject_calc_2(a_f_targ, x_0, Earth, SC)
        optimal_solution_arr.append(solution)
        if solution.craft_obj.t_eject == 1024:
            break
        f = f - f_inc

    print("Done with f_list generation")

    # PLOT??
    t_ej_arr = [solution.craft_obj.t_eject for solution in optimal_solution_arr]
    f_arr = [solution.body_obj.scale_factor for solution in optimal_solution_arr]
    plt.scatter(f_arr, t_ej_arr)
    plt.show()

    start_num = 1
    stop_num = 5
    nums = range(start_num, stop_num+1)
    names = [f'Atm_Data/{num}.txt' for num in nums]
    c1 = 1
    c2 = 1  # Are two values needed? If it's all relative anyway...
    f_dist = np.zeros(len(names))
    # for each pressure profile specified
    for idx, name in enumerate(names):
        print(f'Analyzing Sample: {name}')
        # read in profile, create body
        Earth = copy.copy(Earth_0)
        SC = copy.copy(SC_0)
        Earth.load_density_info(name)  # TODO: fix density func and file loading method
        # vary v, gamma, beta1, beta2
        # find t_ej, a_star(t), a_dot_star(t) (?)
        solution = eject_calc(a_f_targ, x_0, Earth, SC)
        # find distance
        d = np.zeros(len(optimal_solution_arr))
        for idx2 in range(len(optimal_solution_arr)):
            t_step_star = optimal_solution_arr[idx2].t[1] - optimal_solution_arr[idx2].t[0]
            a_star = optimal_solution_arr[idx2].a_t[int(math.floor(optimal_solution_arr[idx2].craft_obj.t_eject / t_step_star))]  # NOTE: relies on constant time steps
            a_dot_star = optimal_solution_arr[idx2].a_dot_t[int(optimal_solution_arr[idx2].craft_obj.t_eject)]
            t_step = solution.t[1] - solution.t[0]
            a = solution.a_t[int(solution.craft_obj.t_eject)]
            a_dot = solution.a_dot_t[int(solution.craft_obj.t_eject)]
            d[idx2] = c1 * (a_star - a)**2.0 + c2 * (a_dot_star - a_dot)**2.0
        # smallest d corresponds to f factor to save
        idx_min = np.argmin(d)
        f_dist[idx] = optimal_solution_arr[idx_min].body_obj.scale_factor

    # show distribution of f factors
    plt.hist(f_dist, len(optimal_solution_arr))  # Could be prettier and centered...
    plt.show()

if __name__ == "__main__":
    main()

    # TODO: must fix density reading and density calculating functions
    # Also needs all around efficiency boosts due to time limits
