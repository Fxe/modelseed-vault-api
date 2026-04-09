"""
Helper class to download data from the Fitness Browser at fit.genomics.lbl.gov.

Usage:
    from fitness_browser import FitnessBrowserOrganism, ORGANISMS

    org = FitnessBrowserOrganism("acidovorax_3H11")
    fitness = org.get_fitness_values()          # raw bytes (tab-delimited)
    t_scores = org.get_t_scores()               # raw bytes (tab-delimited)
    cofit = org.get_cofitness()                  # raw bytes (tab-delimited)
    phenotypes = org.get_specific_phenotypes()   # raw bytes (tab-delimited)
    exp_meta = org.get_experiment_metadata()     # raw bytes (tab-delimited)
    reanno = org.get_reannotations()             # raw bytes (tab-delimited)
    genome = org.get_genome_sequence()           # raw bytes (fasta)
    proteins = org.get_protein_sequences()       # raw bytes (fasta)
    genes = org.get_gene_list()                  # raw bytes (tab-delimited)
"""

import requests
from typing import Optional


BASE_URL = "https://fit.genomics.lbl.gov/cgi-bin"

ORGANISMS = {
    "acidovorax_3H11": "Acidovorax sp. GW101-3H11",
    "azobra": "Azospirillum brasilense Sp245",
    "Btheta": "Bacteroides thetaiotaomicron VPI-5482",
    "Bifido": "Bifidobacterium breve UCC2003",
    "Brev2": "Brevundimonas sp. GW460-12-10-14-LB2",
    "BFirm": "Burkholderia phytofirmans PsJN",
    "Caulo": "Caulobacter crescentus NA1000",
    "Cup4G11": "Cupriavidus basilensis FW507-4G11",
    "PS": "Dechlorosoma suillum PS",
    "DvH": "Desulfovibrio vulgaris Hildenborough JW710",
    "Miya": "Desulfovibrio vulgaris Miyazaki F",
    "Dda3937": "Dickeya dadantii 3937",
    "Ddia6719": "Dickeya dianthicola 67-19",
    "DdiaME23": "Dickeya dianthicola ME23",
    "Dino": "Dinoroseobacter shibae DFL-12",
    "Dyella79": "Dyella japonica UNC79MFTsu3.2",
    "Cola": "Echinicola vietnamensis KMM 6221, DSM 17526",
    "Keio": "Escherichia coli BW25113",
    "HerbieS": "Herbaspirillum seropedicae SmR1",
    "Kang": "Kangiella aquimarina DSM 16071",
    "Koxy": "Klebsiella michiganensis M5al",
    "Lysobacter_OAE881": "Lysobacter sp. OAE881",
    "Magneto": "Magnetospirillum magneticum AMB-1",
    "Marino": "Marinobacter adhaerens HP15",
    "Methanococcus_JJ": "Methanococcus maripaludis JJ",
    "Methanococcus_S2": "Methanococcus maripaludis S2",
    "Mucilaginibacter_YX36": "Mucilaginibacter yixingensis YX-36 DSM 26809",
    "MycoTube": "Mycobacterium tuberculosis H37Rv",
    "Burk376": "Paraburkholderia bryophila 376MFSha3.1",
    "Burkholderia_OAS925": "Paraburkholderia graminis OAS925",
    "Pedo557": "Pedobacter sp. GW460-11-11-14-LB5",
    "Phaeo": "Phaeobacter inhibens DSM 17395",
    "Bvulgatus_CL09T03C04": "Phocaeicola vulgatus CL09T03C04",
    "Ponti": "Pontibacter actiniarum KMM 6156, DSM 19842",
    "pseudo1_N1B4": "Pseudomonas fluorescens FW300-N1B4",
    "pseudo5_N2C3_1": "Pseudomonas fluorescens FW300-N2C3",
    "pseudo6_N2E2": "Pseudomonas fluorescens FW300-N2E2",
    "pseudo3_N2E3": "Pseudomonas fluorescens FW300-N2E3",
    "pseudo13_GW456_L13": "Pseudomonas fluorescens GW456-L13",
    "Putida": "Pseudomonas putida KT2440",
    "WCS417": "Pseudomonas simiae WCS417",
    "psRCH2": "Pseudomonas stutzeri RCH2",
    "SyringaeB728a": "Pseudomonas syringae pv. syringae B728a",
    "SyringaeB728a_mexBdelta": "Pseudomonas syringae pv. syringae B728a \u0394mexB",
    "RalstoniaGMI1000": "Ralstonia solanacearum GMI1000",
    "RalstoniaBSBF1503": "Ralstonia solanacearum IBSBF1503",
    "RalstoniaPSI07": "Ralstonia solanacearum PSI07",
    "RalstoniaUW163": "Ralstonia solanacearum UW163",
    "CL21": "Ralstonia sp. UNC404CL21Col",
    "rhodanobacter_10B01": "Rhodanobacter denitrificans FW104-10B01",
    "RPal_CGA009": "Rhodopseudomonas palustris CGA009",
    "SB2B": "Shewanella amazonensis SB2B",
    "PV4": "Shewanella loihica PV-4",
    "MR1": "Shewanella oneidensis MR-1",
    "ANA3": "Shewanella sp. ANA-3",
    "Smeli": "Sinorhizobium meliloti 1021",
    "Korea": "Sphingomonas koreensis DSMZ 15582",
    "SynE": "Synechococcus elongatus PCC 7942",
    "Variovorax_OAS795": "Variovorax sp. OAS795",
    "Xantho": "Xanthomonas campestris pv. campestris strain 8004",
}


class FitnessBrowserOrganism:
    """Download data from the Fitness Browser for a given organism."""

    def __init__(self, org_id: str, timeout: int = 60):
        if org_id not in ORGANISMS:
            raise ValueError(
                f"Unknown org_id: {org_id!r}. "
                f"Choose from: {', '.join(sorted(ORGANISMS))}"
            )
        self.org_id = org_id
        self.name = ORGANISMS[org_id]
        self.timeout = timeout
        self._session = requests.Session()

    def __repr__(self) -> str:
        return f"FitnessBrowserOrganism({self.org_id!r})  # {self.name}"

    # ── internal ─────────────────────────────────────────────────────

    def _get(self, endpoint: str, params: Optional[dict] = None) -> bytes:
        url = f"{BASE_URL}/{endpoint}"
        p = {"orgId": self.org_id}
        if params:
            p.update(params)
        resp = self._session.get(url, params=p, timeout=self.timeout)
        resp.raise_for_status()
        return resp.content

    # ── public download methods ──────────────────────────────────────

    def get_fitness_values(self) -> bytes:
        """Fitness values (tab-delimited)."""
        return self._get("createFitData.cgi")

    def get_t_scores(self) -> bytes:
        """t scores (tab-delimited)."""
        return self._get("createFitData.cgi", {"t": "1"})

    def get_cofitness(self) -> bytes:
        """Top cofitness for each gene (tab-delimited)."""
        return self._get("createCofitData.cgi")

    def get_specific_phenotypes(self) -> bytes:
        """Specific phenotypes (tab-delimited)."""
        return self._get("spec.cgi", {"download": "1"})

    def get_experiment_metadata(self) -> bytes:
        """Experiment meta-data (tab-delimited)."""
        return self._get("createExpData.cgi")

    def get_reannotations(self) -> bytes:
        """Reannotations with protein sequences (tab-delimited)."""
        return self._get("downloadReanno.cgi")

    def get_genome_sequence(self) -> bytes:
        """Genome sequence (fasta)."""
        return self._get("orgSeqs.cgi", {"type": "nt"})

    def get_protein_sequences(self) -> bytes:
        """Protein sequences (fasta)."""
        return self._get("orgSeqs.cgi")

    def get_gene_list(self) -> bytes:
        """List of genes (tab-delimited)."""
        return self._get("orgGenes.cgi")
