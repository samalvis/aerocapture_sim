# Sam Alvis, January 2022

from common_funcs import *
import copy


def load_lin_atm(filename):
    """Description"""

    return load_outputs(f"Sim_Data/lin_sol_{filename}")


def load_real_atm(filename):
    """Description"""

    return load_outputs(f"Sim_Data/real_sol_{filename}")


def gen_lin_atm(f_inc, f_min, f_max, filename, a_f_targ, body_obj, craft_obj):
    """Description"""

    print("Starting f_list generation")

    lin_solution_arr = []
    f = 1
    # increase f up until crash
    while f <= f_max:
        print(f"Current f value: {f}")
        # adjust body
        body = copy.copy(body_obj)
        craft = copy.copy(craft_obj)
        # find t_ej, a_star(t), a_dot_star(t) (?)
        body.scale_factor = f
        solution = eject_calc(a_f_targ, craft.x_0, body, craft)
        lin_solution_arr.append(solution)
        if solution.craft_obj.t_eject == 0:  # TODO: Fix hardcode? (t_min?)
            break
        f = f + f_inc

    print("Halfway done with f_list generation")

    f = 1 - f_inc
    # decrease f up until escape or 0 atm
    while f >= f_min:
        print(f"Current f value: {f}")
        # adjust body
        body = copy.copy(body_obj)
        craft = copy.copy(craft_obj)
        # find t_ej, a_star(t), a_dot_star(t) (?)
        body.scale_factor = f
        solution = eject_calc(a_f_targ, craft.x_0, body, craft)
        lin_solution_arr.append(solution)
        if solution.craft_obj.t_eject == 1024:  # TODO: Fix hardcode! (t_max?)
            break
        f = f - f_inc

    print("Done with f_list generation")

    # save all linearly varied solutions at once! (need all to work)
    prep_output_file(f"Sim_Data/lin_sol_{filename}")
    save_output(lin_solution_arr, f"Sim_Data/lin_sol_{filename}")

    print("f_list saved for future reference!")

    return lin_solution_arr


def gen_real_atm(start_num, stop_num, filename, body_obj, craft_obj, lin_sol_arr):
    """Description"""

    nums = range(start_num, stop_num+1)
    names = [f'Mars_Atm_Data/{num}.txt' for num in nums]
    c1 = 1
    c2 = 1  # Are two values needed? If it's all relative anyway...
    f_dist = np.zeros(len(names))
    prep_output_file(f"Sim_Data/real_sol_{filename}")  # TODO: make optional? Will overwrite
    # for each pressure profile specified
    for idx, name in enumerate(names):  # Note: could also do multiple reps on each atmosphere, accounting for other variations introduced?
        print(f'Analyzing Sample: {name}')
        # read in profile, create body
        body = copy.deepcopy(body_obj)
        craft = copy.deepcopy(craft_obj)
        body.load_density_info2(name)
        # vary v, gamma, beta1, beta2  # TODO: unhardcode ranges?
        U = craft.x_0[3]  # TODO: must vary V? not U?
        craft.x_0[3] = np.random.normal(U, 0.01/3.0)
        gamma = craft.x_0[4]
        craft.x_0[4] = np.random.normal(gamma, math.radians(0.2/3.0))
        b_i = craft.beta_i
        craft.beta_i = np.random.uniform(0.95 * b_i, 1.05 * b_i)
        b_f = craft.beta_f
        craft.beta_f = np.random.uniform(0.95 * b_f, 1.05 * b_f)
        # find t_ej for this case
        solution = eject_calc(a_f_targ, craft.x_0, body, craft)
        plot_case(solution)
        a = 0
        a_dot = 0
        if solution.craft_obj.t_eject > solution.t[-1]:
            print("Impact?")
            print(solution.craft_obj.x_0)
            print(solution.craft_obj.t_eject)
            print(solution.t[-1])
            a = solution.a_t[-1]
            a_dot = solution.a_dot_t[-1]
        else:
            a = solution.a_t[int(solution.craft_obj.t_eject)]
            a_dot = solution.a_dot_t[int(solution.craft_obj.t_eject)]

        # find distance to linear atmosphere t_ej match candidates
        d = np.zeros(len(lin_sol_arr))
        for idx2 in range(len(lin_sol_arr)):
            t_step_star = lin_sol_arr[idx2].t[1] - lin_sol_arr[idx2].t[0]
            a_star = lin_sol_arr[idx2].a_t[int(math.floor(lin_sol_arr[idx2].craft_obj.t_eject / t_step_star))]  # NOTE: relies on constant time steps
            a_dot_star = lin_sol_arr[idx2].a_dot_t[int(lin_sol_arr[idx2].craft_obj.t_eject)]
            d[idx2] = c1 * (a_star - a)**2.0 + c2 * (a_dot_star - a_dot)**2.0
        # smallest d corresponds to f factor to save
        idx_min = np.argmin(d)
        print(lin_sol_arr[idx_min].body_obj.scale_factor)
        f_dist[idx] = lin_sol_arr[idx_min].body_obj.scale_factor
        save_output(f_dist[idx], f"Sim_Data/real_sol_{filename}")
        # save real solution to solutions array?

    return f_dist


def plot_f_dist(f_dist, n_bins):  # n_bins from number of
    """Description"""

    plt.hist(f_dist, bins=n_bins)  # Could be prettier and centered...
    plt.show()


if __name__ == "__main__":
    # define body
    Mars = Body()
    Mars.mu = 4.2828E4
    Mars.R = 3396.2
    Mars.J2 = 1.9555E-3
    Mars.omega = 2 * math.pi / 88619.6
    Mars.max_atm_height = 125.0
    Mars.load_density_info2("Mars_Atm_Data/0.txt")

    # define craft
    SC = Craft()
    SC.beta_i = 30.0
    SC.beta_f = 135.0
    SC.t_eject = 0  # NOTE: here, this time will be adjusted and so is irrelevant

    # define target apoapsis
    a_f_targ = 500 + Mars.R

    # define x_0
    SC.x_0 = [Mars.R + 125.0, math.radians(18.44), math.radians(77.45), 6, 0, math.radians(45)]
    # SC.x_0[4] = EFPA_calc(a_f_targ, x_0, Mars, SC).craft_obj.x_0[4]
    SC.x_0[4] = math.radians(-9.875)

    # load linear distributions
    lin_sols = load_lin_atm("test5")

    # OR generate distributions
    # lin_sols = gen_lin_atm(0.01, 0.75, 1.25, "test4", a_f_targ, Mars, SC)

    # load atmospheres
    # f_dist = load_real_atm("test8_var")

    # OR generate atmospheres (remember to disperse!)
    f_dist = gen_real_atm(1, 20, "test10_var", Mars, SC, lin_sols)  # return solutions? which?

    # plot results
    plot_f_dist(f_dist, len(lin_sols))
