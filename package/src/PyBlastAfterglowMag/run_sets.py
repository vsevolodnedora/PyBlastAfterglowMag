import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from .interface import PyBlastAfterglow, Ejecta

class CollateDataForRuns:

    def __init__(self, workingdirs : list, key_to_key_map : dict, verbose : bool = False):

        self.workingdirs = workingdirs
        if len(self.workingdirs) == 0:
            raise ValueError("No workingdirs are given")

        self.key_to_key_map = key_to_key_map
        if len(key_to_key_map.keys()) == 0:
            raise ValueError("Empty 'from_to_keys' dictionary is given")

        self.out = {}
        for in_key, out_key in self.key_to_key_map.items():
            self.out[out_key] = []

        self.verb = verbose

    def _process_workdir(self, i : int, KL : Ejecta, do_id : bool, do_lc : bool, do_skymap : bool) -> None:
        """
        Load files for a single simulation and add the necessary quantities to the dataframe
        :param i:
        :param KL:
        :param do_id:
        :param do_lc:
        :param do_skymap:
        :return:
        """

        if do_skymap:
            raise NotImplementedError("not implemented yet")

        if do_lc:
            times = KL.get_lc_times(spec=False)
            freqs = KL.get_lc_freqs(spec=False, unique=True)
            lc_attrs = KL.get_lc_obj().attrs
            for freq in freqs:
                for time in times:

                    for in_key, out_key in self.key_to_key_map.items():
                        if (in_key in lc_attrs.keys()):
                            self.out[out_key].append(np.float64(lc_attrs[in_key]))

                    if do_id:
                        id_attrs = KL.get_id_obj().attrs
                        for in_key, out_key in self.key_to_key_map.items():
                            if ((in_key in id_attrs.keys()) and (not in_key in lc_attrs.keys())):
                                self.out[out_key].append( id_attrs[in_key] )


                    flux = KL.get_lc_totalflux(freq=freq,time=time,spec=False)
                    if "time" in self.key_to_key_map.keys():
                        self.out[self.key_to_key_map["time"]].append( np.float64(time) )
                    if "freq" in self.key_to_key_map.keys():
                        self.out[self.key_to_key_map["freq"]].append( np.float64(freq) )
                    if "flux" in self.key_to_key_map.keys():
                        self.out[self.key_to_key_map["flux"]].append( np.float64(flux) )

            for _, key in self.key_to_key_map.items():
                if (len(self.out[key])) != len(times) * len(freqs) * (i + 1):
                    raise ValueError("Size mismatch")

    def process(self, mode : str, do_id : bool, do_lc : bool, do_skymap : bool):
        for i in tqdm(range(len(self.workingdirs))):
            if (not os.path.isdir(self.workingdirs[i])):
                raise IOError(f"Directory is not found: {self.workingdirs[i]}")

            # initialize PBA instance in the working directory
            workingdir = self.workingdirs[i]
            pba = PyBlastAfterglow(workingdir=workingdir,verbose=self.verb)
            if mode == "GRB": KL = pba.GRB
            elif mode == "KN": KL = pba.KN
            else: raise NotImplementedError(f"mode = {mode} is not supported")

            self._process_workdir(i=i, KL=KL, do_id=do_id, do_lc=do_lc, do_skymap=do_skymap)

        print(f"Processing {len(self.workingdirs)} complete")

    def get_df(self) -> pd.DataFrame:
        df = pd.DataFrame.from_dict(self.out)
        return df


