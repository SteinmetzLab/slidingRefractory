from pathlib import Path
import numpy as np
import pandas as pd

class SteinmetzDataAccess:
    
    def __init__(self, clusters_path, data_dir):
        self.clusters_path = clusters_path
        self.data_dir = data_dir
        self._load()
        
    def _load(self):
        self.clusters = pd.read_parquet(self.clusters_path)
        self.insertions = self.clusters.index.get_level_values(0).unique()
        
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
            cluster_idx_reg = ins_table.index.to_numpy()
            acgs, cluster_ids = self.load_acgs_insertion(ins, kind=kind)
            # ensure we're grabbing the correct units
            assert np.all(cluster_ids[cluster_idx_reg] == ins_table.cluster_id.to_numpy())
            acgs = acgs[:, cluster_idx_reg]
            acg_list.append(acgs)
        out_acgs = np.concatenate(acg_list, axis=1)
        
        assert out_acgs.shape[1] == len(table)
        
        return out_acgs, table
            
        
    def list_subjects(self):
        """
        List Steinmetz subjects.
        """
        return list(np.unique([self._split_steinmetz_insertion(st)[0] for st in self.insertions]))
    
    def list_sessions(self, subject):
        """
        List sessions for a given subject.
        
        :param subject: Subject name
        """
        subj = [p for p in self.insertions if p.startswith(subject)]
        ses = np.unique([self._split_steinmetz_insertion(st)[1] for st in subj])
        return list(ses)
    
    def list_probes(self, subject, session):
        """
        List probes for a given subject and session.
        
        :param subject: Subject name
        :param session: Session (date)
        """
        subj = [p for p in self.insertions if p.startswith(subject)]
        ses = [p for p in subj if p.split("_")[1] == session]
        prb = np.unique([self._split_steinmetz_insertion(st)[2] for st in ses])
        return list(prb)
    
    def list_insertions(self):
        """
        List all insertions.
        """
        return self.insertions
    
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
            
        subject, session, probe = self._split_steinmetz_insertion(insertion)
        ins_path = self.data_dir.joinpath(subject, session, "ephys_" + probe)
        acg = np.load(ins_path.joinpath(f"correlograms_{string}.npy"))
        cluster_ids = np.load(ins_path.joinpath("cluster_ids.npy"))
        
        return acg, cluster_ids
        
    def _split_steinmetz_insertion(self, ins):
        """
        Splits insertion name for path finding.
        """
        split = ins.split("_")
        subject = split[0]
        session = split[1]
        probe = split[2]
        return subject, session, probe