"""
Helper class to download data from the Fitness Browser at fit.genomics.lbl.gov.

The site is behind Cloudflare Bot Management. This module keeps a real Chrome
browser open (via undetected-chromedriver) and performs all downloads through it.

Setup:
    pip install undetected-chromedriver

Usage:
    from fitness_browser import FitnessBrowser, ORGANISMS

    with FitnessBrowser(headless=False, version_main=146) as fb:
        # Download data for any organism
        fitness = fb.get_fitness_values("acidovorax_3H11")
        genes = fb.get_gene_list("Keio")
        proteins = fb.get_protein_sequences("MR1")

        # Loop over all organisms
        for org_id in ORGANISMS:
            data = fb.get_fitness_values(org_id)
"""

import logging
import time
from typing import Optional
from urllib.parse import quote

logger = logging.getLogger(__name__)

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


def _wait_for_cloudflare(driver, timeout: int = 30):
    """Wait for Cloudflare challenge to resolve."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        title = driver.title or ""
        if "Just a moment" not in title and "Attention" not in title:
            if title:
                return
        time.sleep(1)
    raise TimeoutError(
        f"Cloudflare challenge not solved within {timeout}s. "
        f"Try headless=False to see what's happening."
    )


class FitnessBrowser:
    """Single browser session for downloading from the Fitness Browser.

    Opens Chrome once, solves Cloudflare once, then downloads data for
    any organism through the same session.

    Args:
        headless: Run Chrome headless (no visible window). Set False to debug.
        version_main: Force Chrome major version for ChromeDriver (e.g. 146).
                      Find yours at chrome://version.
        timeout: Timeout in seconds for Cloudflare challenge and downloads.
    """

    def __init__(
        self,
        headless: bool = True,
        version_main: Optional[int] = None,
        timeout: int = 60,
    ):
        self.timeout = timeout
        self._driver = None
        self._start_browser(headless=headless, version_main=version_main)

    def _start_browser(self, headless: bool, version_main: Optional[int]):
        """Launch Chrome and solve the Cloudflare challenge."""
        import undetected_chromedriver as uc

        options = uc.ChromeOptions()
        if headless:
            options.add_argument("--headless=new")
        options.add_argument("--no-sandbox")
        options.add_argument("--disable-dev-shm-usage")

        self._driver = uc.Chrome(options=options, version_main=version_main)

        # Navigate to the main page to solve Cloudflare
        logger.info("Solving Cloudflare challenge...")
        self._driver.get(f"{BASE_URL}/orgAll.cgi")
        _wait_for_cloudflare(self._driver, timeout=self.timeout)
        logger.info("Cloudflare challenge solved.")

    def __repr__(self) -> str:
        return f"FitnessBrowser(driver={'open' if self._driver else 'closed'})"

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        """Close the browser."""
        if self._driver:
            try:
                self._driver.quit()
            except Exception:
                pass
            self._driver = None

    # -- internal --

    def _validate_org_id(self, org_id: str):
        if org_id not in ORGANISMS:
            raise ValueError(
                f"Unknown org_id: {org_id!r}. "
                f"Choose from: {', '.join(sorted(ORGANISMS))}"
            )

    def _fetch(self, url: str) -> bytes:
        """Fetch a URL using the browser's own fetch() API."""
        script = """
        async function fetchData(url) {
            const resp = await fetch(url, {credentials: 'include'});
            if (!resp.ok) {
                throw new Error(`HTTP ${resp.status}: ${resp.statusText}`);
            }
            const buffer = await resp.arrayBuffer();
            const bytes = new Uint8Array(buffer);
            let binary = '';
            const chunkSize = 8192;
            for (let i = 0; i < bytes.length; i += chunkSize) {
                const chunk = bytes.subarray(i, i + chunkSize);
                binary += String.fromCharCode.apply(null, chunk);
            }
            return btoa(binary);
        }
        return await fetchData(arguments[0]);
        """
        import base64

        try:
            b64_data = self._driver.execute_script(script, url)
        except Exception as e:
            raise RuntimeError(f"Browser fetch failed for {url}: {e}") from e

        return base64.b64decode(b64_data)

    def _get(self, org_id: str, endpoint: str, params: Optional[dict] = None) -> bytes:
        self._validate_org_id(org_id)
        p = {"orgId": org_id}
        if params:
            p.update(params)
        query = "&".join(f"{k}={quote(str(v))}" for k, v in p.items())
        url = f"{BASE_URL}/{endpoint}?{query}"
        return self._fetch(url)

    # -- public download methods --

    def get_fitness_values(self, org_id: str) -> bytes:
        """Fitness values (tab-delimited)."""
        return self._get(org_id, "createFitData.cgi")

    def get_t_scores(self, org_id: str) -> bytes:
        """t scores (tab-delimited)."""
        return self._get(org_id, "createFitData.cgi", {"t": "1"})

    def get_cofitness(self, org_id: str) -> bytes:
        """Top cofitness for each gene (tab-delimited)."""
        return self._get(org_id, "createCofitData.cgi")

    def get_specific_phenotypes(self, org_id: str) -> bytes:
        """Specific phenotypes (tab-delimited)."""
        return self._get(org_id, "spec.cgi", {"download": "1"})

    def get_experiment_metadata(self, org_id: str) -> bytes:
        """Experiment meta-data (tab-delimited)."""
        return self._get(org_id, "createExpData.cgi")

    def get_reannotations(self, org_id: str) -> bytes:
        """Reannotations with protein sequences (tab-delimited)."""
        return self._get(org_id, "downloadReanno.cgi")

    def get_genome_sequence(self, org_id: str) -> bytes:
        """Genome sequence (fasta)."""
        return self._get(org_id, "orgSeqs.cgi", {"type": "nt"})

    def get_protein_sequences(self, org_id: str) -> bytes:
        """Protein sequences (fasta)."""
        return self._get(org_id, "orgSeqs.cgi")

    def get_gene_list(self, org_id: str) -> bytes:
        """List of genes (tab-delimited)."""
        return self._get(org_id, "orgGenes.cgi")
