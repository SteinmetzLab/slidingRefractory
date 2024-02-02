from pathlib import Path
import numpy as np
import pandas as pd

class AllenDataAccess:
    
    def __init__(self, clusters_path, data_dir):
        self.clusters_path = clusters_path
        self.data_dir = data_dir
        self._load()
        
    def _load(self):
        self.clusters = pd.read_parquet(self.clusters_path)
        self.insertions = self.clusters.index.get_level_values(0).unique()
        
    def list_insertions(self):
        """
        List all insertions.
        """
        return self.insertions
    
    def session2probes(self, session):
        return list(clusters[clusters.ecephys_session_id==sid].index.get_level_values(0).unique())
        
    def probe2session(self, probe):
        return self.clusters.loc[probe].ecephys_session_id.iloc[0]
        
    def load_acgs_region(self, region, kind="rf"):
        """
        Load ACGs and metrics for units from a given cosmos acronym.
        
        :param region: Cosmos acronym
        :param kind: "rf" or "ts"
        :return: acgs (np array), table (DataFrame)
        """
        assert kind in ["rf", "ts"], "kind must be rf or ts."
        if kind == "rf":
            string = "refractory_period"
        elif kind == "ts":
            string = "time_scale"
            
        table = self.clusters[self.clusters.cosmos_acronym == region]
        acg_list = []
        for ins in table.index.get_level_values(0).unique():
            ins_table = table.loc[ins]
            cluster_id_reg = ins_table.index.to_numpy()
            acgs, cluster_ids = self.load_acgs_insertion(ins, kind=kind)
            cluster_idx_reg = np.searchsorted(cluster_ids, cluster_id_reg)
            # ensure we're grabbing the correct units
            assert np.all(cluster_ids[cluster_idx_reg] == ins_table.index.to_numpy())
            acgs = acgs[:, cluster_idx_reg]
            acg_list.append(acgs)
        out_acgs = np.concatenate(acg_list, axis=1)
        
        assert out_acgs.shape[1] == len(table)
        
        return out_acgs, table
    
    def load_acgs_insertion(self, insertion, kind="rf"):
        """
        Returns acg array for this insertion along with cluster ids.
        :param insertion: Name of insertion
        :param kind: 'rf' or 'ts'
        :return: acg, cluster_ids
        """
        assert kind in ["rf", "ts"], "kind must be rf or ts."
        if kind == "rf":
            string = "refractory_period"
        elif kind == "ts":
            string = "time_scale"
        
        session_id = self.probe2session(insertion)
        ins_path = self.data_dir.joinpath(str(session_id), str(insertion))
        acg = np.load(ins_path.joinpath(f"correlograms_{string}.npy"))
        cluster_ids = np.load(ins_path.joinpath("cluster_ids.npy"))
        
        return acg, cluster_ids