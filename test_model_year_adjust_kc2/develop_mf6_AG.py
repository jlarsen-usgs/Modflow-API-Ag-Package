from modflowapi import ModflowApi
import numpy as np
import pandas as pd
import flopy
import os


def create_test_model(name):
    import flopy as fp

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



class Modflow6Ag(object):
    """

    """
    def __init__(self, sim):

        self.sim = sim
        name = list(sim.model_names)[0]
        self.gwf = sim.get_model(name)
        self.name = name.upper()
        self.ncpl = self.gwf.modelgrid.ncpl
        self.nlay = self.gwf.modelgrid.nlay
        self.totim = np.add.accumulate(self.gwf.modeltime.perlen)
        self.area = None
        self.top = None
        self.top_tile = None
        self.botm = None

    def create_addresses(self, mf6):
        sto_name = self.gwf.sto.name[0].upper()
        dis_name = self.gwf.dis.name[0].upper()
        try:
            rch_name = self.gwf.rcha.name[0].upper()
        except AttributeError:
            rch_name = self.gwf.rch.name[0].upper()

        try:
            evt_name = self.gwf.evta.name[0].upper()
        except AttributeError:
            evt_name = self.gwf.evta.name[0].upper()

        self.area = mf6.get_value(
            mf6.get_var_address("AREA", self.name, dis_name)
        )

        self.sy = mf6.get_value(
            mf6.get_var_address("SY", self.name, sto_name)
        )
        self.ss = mf6.get_value(
            mf6.get_var_address("SS", self.name, sto_name)
        )
        self.botm = mf6.get_value(
            mf6.get_var_address("BOT", self.name, dis_name)
        )

        self.top = mf6.get_value(
            mf6.get_var_address("TOP", self.name, dis_name)
        )

        self.top_tile = np.tile(self.top, self.nlay)

        # todo: set the well name!!!!
        self.iusesy_addr = mf6.get_var_address("IUSESY", self.name, sto_name)
        self.well_addr = mf6.get_var_address("BOUND", self.name, "WEL_0")
        self.well_node_addr = mf6.get_var_address(
            "NODELIST", self.name, "WEL_0"
        )
        self.rch_addr = mf6.get_var_address("BOUND", self.name, rch_name)
        self.evt_addr = mf6.get_var_address("BOUND", self.name, evt_name)
        self.head_addr = mf6.get_var_address("X", self.name)

    def calculate_demand(self, mf6, well_ids, kc):
        # todo: change well_ids to crop_nodes!
        well_nodes = mf6.get_value(self.well_node_addr)
        et = mf6.get_value(self.evt_addr)

        node = well_nodes[well_ids]
        demand = (et[node, 1] * kc)

        # adjust the et value so we don't double account
        # todo: actually I don't think that this is quite correct...
        #  we should keep et as is and then adjust the demand value.
        #  more or less provide accounting of the ET, and then let
        #  other modflow processes deal with it!
        if demand > et[node, 1]:
            et[node, 1] = 0
        else:
            et[node, 1] -= demand

        mf6.set_value(self.evt_addr, et)
        return demand

    def calculate_gw_demand(self, mf6, demand, crop_node, rooting_depth):
        # todo: need rooting depth for crop nodes
        # will be something like if head > layer top - rooting depth
        # else gw_volume = 0

        head = mf6.get_value(self.head_addr)
        iusesy = mf6.get_value(self.iusesy_addr)



        # gw_volume needs to be changed to a rooting depth calc.
        gw_column = (head - (self.top - rooting_depth))
        gw_volume = gw_column * self.area
        # gw_volume = (head - self.botm) * self.area
        gw_volume[gw_volume < 0] = 0
        gw_volume = np.where(
            iusesy > 0,
            gw_volume * self.sy,
            0 # gw_volume * self.ss
        )

        # todo: calculate the drop in elevation not the absoulte elevation.
        q = gw_volume[crop_node] - demand
        if q > 0:
            # take supply from groundwater!
            gw_volume[crop_node] = q
            gw_demand = demand
        else:
            gw_demand = gw_volume[crop_node]
            gw_volume[crop_node] = 0

        gw_head = np.where(
            iusesy > 0,
            gw_volume / self.sy,
            gw_volume / self.ss
        )

        adj_gw_column = np.where(
            iusesy > 0,
            gw_volume / self.sy,
            gw_volume / self.ss
        ) / self.area

        print(gw_column, adj_gw_column)

        head_adj = np.where(adj_gw_column > 0,
                            gw_column - adj_gw_column,
                            0)

        ohead = head - head_adj
        # ohead = np.where(gw_volume < 0,
        #                  self.botm,
        #                  gw_head / self.area)
        print(f"{head[0] :.3f} > {ohead[0] :.3f}")

        mf6.set_value(self.head_addr, ohead)

        return gw_demand

    def calculate_unsaturated_demand(self, mf6, demand, gw_demand, crop_node):
        # todo: need to implement...

        return 0

    def calculate_precipiation_demand(
        self, mf6, demand, gw_demand, uzf_demand, crop_node
    ):
        ppt = mf6.get_value(self.rch_addr)

        q = (ppt[crop_node] + uzf_demand + gw_demand) - demand

        if q > 0:
            ppt_demand = demand - (uzf_demand + gw_demand)
            ppt[crop_node] = q
        else:
            ppt_demand = ppt[crop_node]
            ppt[crop_node] = 0

        mf6.set_value(self.rch_addr, ppt)
        return ppt_demand

    def calculate_pumping_demand(
            self, mf6, demand, gw_demand, uzf_demand, ppt_demand, well_id
    ):
        # well_nodes = mf6.get_value(self.well_node_addr)
        # node = well_nodes[well_id]
        well = mf6.get_value(self.well_addr)
        q = demand - (gw_demand + uzf_demand + ppt_demand)

        if q > np.abs(well[0, well_id]):
            pumping_demand = well[0, well_id]
        else:
            well[0, well_id] = -q
            pumping_demand = q

        mf6.set_value(self.well_addr, well)

        return pumping_demand

    def run_model(self, dll):

        mf6 = ModflowApi(dll)
        mf6.initialize()
        current_time = mf6.get_current_time()
        end_time = mf6.get_end_time()

        self.create_addresses(mf6)
        # todo: add this to the create_addresses routine
        max_iter = mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))

        # prepare the iteration loops
        kper = 0
        actual_demand = []
        while current_time < end_time:
            dt = mf6.get_time_step()
            mf6.prepare_time_step(dt)
            kiter = 0

            if current_time in self.totim or current_time == 0.:
                if current_time in self.totim:
                    kper += 1

                # todo: update this for some sort of input (crop id, kc value)
                #  wrt time stepping: demand is not going to change, but
                #  gw_demand, uzf, ppt, and pumping demand should be updated
                #  each ts.
                demand = self.calculate_demand(mf6, 0, 1.0)

                # todo: update this for some sort of crop id input
                #  and rooting depth array, this should most likely be calculated
                #  each and every time step!
            gw_demand = self.calculate_gw_demand(mf6, demand, 0, 6)
            uzf_demand = self.calculate_unsaturated_demand(
                mf6, demand, gw_demand, 0
            )
            ppt_demand = self.calculate_precipiation_demand(
                mf6, demand, gw_demand, uzf_demand, 0
            )
            # todo: make sure that this gets a well node!!!
            pumping_demand = self.calculate_pumping_demand(
                mf6, demand, gw_demand, uzf_demand, ppt_demand, 0
            )

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

        return success

if __name__ == "__main__":
    dll = os.path.join("..", "modflow-bmi", "libmf6.dll")
    sim, gwf = create_test_model("GWF")

    mf6ag = Modflow6Ag(sim)
    mf6ag.run_model(dll, )

