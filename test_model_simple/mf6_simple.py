from modflowapi import ModflowApi
import os


def create_test_model(name):
    import flopy as fp

    sim_ws = os.path.join(".",)
    sim = fp.mf6.MFSimulation(name, sim_ws=sim_ws)

    # create TDIS
    period_data = [(i + 1, 1, 1.0) for i in range(10)]
    tdis = fp.mf6.ModflowTdis(sim, nper=10, perioddata=tuple(period_data))

    # create IMS
    ims = fp.mf6.ModflowIms(sim)

    # create model!
    gwf = fp.mf6.ModflowGwf(sim, modelname="simple_ex", save_flows=True, print_input=True, print_flows=True)

    dis = fp.mf6.ModflowGwfdis(
        gwf,
        nrow=10,
        ncol=10,
        delr=10,
        delc=10,
        top=100
    )

    ic = fp.mf6.ModflowGwfic(gwf, strt=100)
    npf = fp.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True)

    stress_period_data = {
        i: [[(0, 4, 4), 0.2], [(0, 9, 9), 0.2]] for i in range(10)
    }
    wel = fp.mf6.ModflowGwfwel(gwf, stress_period_data=stress_period_data)

    budget_file = f"{name}.cbc"
    head_file = f"{name}.hds"
    saverecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(10)}
    printrecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(10)}
    oc = fp.mf6.ModflowGwfoc(gwf,
                             budget_filerecord=budget_file,
                             head_filerecord=head_file,
                             saverecord=saverecord,
                             printrecord=printrecord)

    sim.write_simulation()

    return sim, gwf


if __name__ == "__main__":
    gwfname = "GWF"
    dll = os.path.join("..", "..", "modflow-bmi", "libmf6.dll")

    sim, gwf = create_test_model(gwfname)

    mf6 = ModflowApi(dll)
    mf6.initialize()

    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    max_iter = mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))

    well_addr = mf6.get_var_address("BOUND", gwfname.upper(), "WEL")

    # prepare the iteration loops
    while current_time < end_time:
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)
        kiter = 0
        well_flux = mf6.get_value(well_addr)
        n_solutions = mf6.get_subcomponent_count()
        for sol_id in range(1, n_solutions + 1):

            mf6.prepare_solve(sol_id)
            while kiter < max_iter:
                has_converged = mf6.solve(sol_id)
                kiter += 1
                if has_converged:
                    break

            mf6.finalize_solve(sol_id)

        mf6.finalize_time_step()
        current_time = mf6.get_current_time()

        if not has_converged:
            print("model did not converge")

    try:
        mf6.finalize()
        success = True
    except:
        raise RuntimeError()

    budget = gwf.output.budget()
    rec_names = budget.get_unique_record_names()
    wel = budget.get_data(text=rec_names[-1])
    print(wel)
