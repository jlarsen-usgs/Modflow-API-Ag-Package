from modflowapi import ModflowApi
import os


def create_test_model(name):
    import flopy as fp
    import pandas as pd

    sim_ws = os.path.join(".",)
    sim = fp.mf6.MFSimulation(name, sim_ws=sim_ws)

    # create TDIS with monthly stress periods and daily time steps
    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    period_data = [(i, i, 1.0) for i in perlen]
    tdis = fp.mf6.ModflowTdis(
        sim,
        nper=12,
        perioddata=tuple(period_data),
        time_units="days"
    )

    # create IMS
    ims = fp.mf6.ModflowIms(sim, complexity="MODERATE")

    # create model!
    gwf = fp.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        print_input=True,
        print_flows=True
    )

    # define delc and delr to equal approximately 1 acre
    dis = fp.mf6.ModflowGwfdis(
        gwf,
        nrow=10,
        ncol=10,
        delr=63.6,
        delc=63.6,
        top=100,
        length_units='meters'
    )

    ic = fp.mf6.ModflowGwfic(gwf, strt=95)
    npf = fp.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True)
    sto = fp.mf6.ModflowGwfsto(gwf, iconvert=1)

    stress_period_data = {
        i: [[(0, 4, 4), -100.], [(0, 9, 9), -100.]] for i in range(12)
    }
    wel = fp.mf6.ModflowGwfwel(gwf, stress_period_data=stress_period_data)

    # create RCH and EVT packages from DAVIS monthly average CIMIS data
    cimis_data = os.path.join("..", "data", "davis_monthly_ppt_eto.txt")
    df = pd.read_csv(cimis_data)

    recharge = {i: v / perlen[i] for i, v in enumerate(df.ppt_avg_m.values)}
    rch = fp.mf6.ModflowGwfrcha(gwf, recharge=recharge)

    surface = {i: 100 for i in range(12)}
    eto = {i: v / perlen[i] for i, v in enumerate(df.eto_avg_m.values)}
    depth = {i: 3 for i in range(12)}
    evt = fp.mf6.ModflowGwfevta(gwf, surface=surface, rate=eto, depth=depth)

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
    dll = os.path.join("..", "modflow-bmi", "libmf6.dll")

    sim, gwf = create_test_model(gwfname)

    mf6 = ModflowApi(dll)
    mf6.initialize()

    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    max_iter = mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))

    well_addr = mf6.get_var_address("BOUND", gwfname.upper(), "WEL_0")
    well_node_addr = mf6.get_var_address("NODELIST", gwfname.upper(), "WEL_0")
    rch_addr = mf6.get_var_address("BOUND", gwfname.upper(), "RCHA_0")
    evt_addr = mf6.get_var_address("BOUND", gwfname.upper(), "EVTA_0")

    # prepare the iteration loops
    while current_time < end_time:
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)
        kiter = 0

        well_value = mf6.get_value(well_addr)
        well_nodes = mf6.get_value(well_node_addr)
        rch_value = mf6.get_value(rch_addr)
        evt_value = mf6.get_value(evt_addr)

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
    wel = budget.get_data(text="WEL")
    print(rec_names)
    print(wel)
