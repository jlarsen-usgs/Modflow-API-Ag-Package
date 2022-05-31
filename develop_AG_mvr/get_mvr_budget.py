import numpy as np
import pandas as pd
import flopy
import os


headings = {
    "time_step": -999,
    "stress_period": -999,
    "totim": -999,
    "tslen": -999,
    "provider": 0,
    "pid": 1,
    "qa": 2,
    "qp": 3,
    "reciever": 4,
    "rid": 5
}

class MvrBudget(flopy.utils.mflistfile.ListBudget):
    """

    """
    def __init__(self, file_name, budgetkey=None, timeunit="days"):

        super().__init__(file_name, budgetkey, timeunit)

    def set_budget_key(self):
        self.budgetkey = "WATER MOVER PACKAGE (MVR) FLOW RATES"
        return

    def _load(self, maxentries=None):
        self._build_index(maxentries)
        totim = []
        i0 = 0
        i1 = 0
        d = {i: [] for i in headings.keys()}
        for ts, sp, seekpoint in self.idx_map:
            tinc, tcum = self._get_sp(ts, sp, seekpoint)
            for rec in tinc:
                for key, index in headings.items():
                    if index == -999:
                        continue
                    d[key].append(rec[index])
                i0 += 1

            # Get the time for this record
            seekpoint = self._seek_to_string("TIME SUMMARY AT END")
            tslen, sptim, tt = self._get_totim(ts, sp, seekpoint)
            totim.append(tt)

            while i1 < i0:
                d["time_step"].append(ts)
                d["stress_period"].append(sp)
                d["totim"].append(tt)
                d["tslen"].append(tslen)
                i1 += 1

        inc = pd.DataFrame.from_dict(d)


        for col in list(inc):
            try:
                inc[col] = pd.to_numeric(inc[col])
            except ValueError:
                pass

        inc["time_step"] -= 1
        inc["stress_period"] -= 1
        inc["pid"] -= 1
        inc["rid"] -= 1

        self.inc = inc
        self.cum = inc.copy()

    def _get_index(self, maxentries):
        # --parse through the file looking for matches and parsing ts and sp
        idxs = []
        l_count = 1
        while True:
            seekpoint = self.f.tell()
            line = self.f.readline()
            if line == "":
                break
            if self.budgetkey in line:
                for _ in range(self.tssp_lines):
                    line = self.f.readline()
                try:
                    ts, sp = get_ts_sp(line)
                except:
                    print(
                        "unable to cast ts,sp on line number",
                        l_count,
                        " line: ",
                        line,
                    )
                    break
                # print('info found for timestep stress period',ts,sp)

                idxs.append([ts, sp, seekpoint])

                if maxentries and len(idxs) >= maxentries:
                    break

        return idxs

    def _get_sp(self, ts, sp, seekpoint):
        self.f.seek(seekpoint)
        n = 0
        recs = []
        while True:
            line = self.f.readline()
            if n < 5:
                n += 1
                continue
            if line == "":
                print(
                    "end of file found while seeking budget "
                    "information for ts,sp: {} {}".format(ts, sp)
                )
                return self.null_entries

            elif "----------------" in line:
                break

            rec = line.strip().split()[1:]
            recs.append(rec)

        return recs, recs


def get_ts_sp(line):
    line = line.replace(",", "").replace("*", "")

    searchstring = "STEP"
    idx = line.index(searchstring) + len(searchstring)
    ll = flopy.utils.flopy_io.line_parse(line[idx:])
    ts = int(ll[0])

    searchstring = "PERIOD"
    idx = line.index(searchstring) + len(searchstring)
    ll = flopy.utils.flopy_io.line_parse(line[idx:])
    sp = int(ll[0])

    return ts, sp
