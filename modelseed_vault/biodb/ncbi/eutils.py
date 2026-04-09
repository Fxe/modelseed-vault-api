import logging
import os
import requests
import xml.etree.ElementTree as ET


logger = logging.getLogger(__name__)


def parse_link_set_db(tree):
    res = set()
    for el in tree:
        if el.tag == "Link":
            for o in el:
                if o.tag == "Id":
                    res.add(int(o.text))
    return res


def parse_link_set(tree):
    link_set = {}
    for el in tree:
        if el.tag == "LinkSetDb":
            link_set["link_set_db"] = parse_link_set_db(el)
    return link_set


def get_link_set(res):
    link_set = None
    tree = ET.fromstring(res.content)
    for el in tree:
        if el.tag == "LinkSet":
            link_set = parse_link_set(el)
    return link_set


def parse_document_summary_set(tree):
    documents = []
    for el in tree:
        if el.tag == "DocumentSummary":
            documents.append(parse_document_summary(el))
        else:
            print("DocumentSummary?", el, el.tag, el.text, el.attrib)
    return documents


def _to_item(e):
    doc_item = dict(e.items())
    doc_item['_tag'] = e.tag
    if e.text:
        doc_item['_text'] = e.text
    return doc_item


def parse_biosample(t):
    """
    tree = ET.fromstring(xml_str)
    docs = []
    for el in tree:
        if el.tag == "BioSample":
            docs.append(parse_biosample(el))
        else:
            print(el.tag, el)
    :param t:
    :return:
    """
    doc_sum = {}
    for el in t:
        if el.tag == "Ids":
            doc_sum["Ids"] = []
            for e in el:
                doc_sum["Ids"].append(_to_item(e))
        elif el.tag == "Description":
            doc_sum["Description"] = []
            for e in el:
                if e.tag == 'Organism':
                    org = _to_item(e)
                    for _e in e:
                        org[f'_{_e.tag}'] = _to_item(_e)
                    doc_sum["Description"].append(org)
                else:
                    doc_sum["Description"].append(_to_item(e))
        elif el.tag == "Owner":
            doc_sum["Owner"] = []
            for e in el:
                if e.tag == 'Contacts':
                    ct = _to_item(e)
                    for _e in e:
                        ct[f'_{_e.tag}'] = _to_item(_e)
                        for __e in _e:
                            ct[f'_{_e.tag}'][f'_{__e.tag}'] = _to_item(__e)
                            for ___e in __e:
                                ct[f'_{_e.tag}'][f'_{__e.tag}'][f'_{___e.tag}'] = _to_item(___e)
                    doc_sum["Owner"].append(ct)
                else:
                    doc_sum["Owner"].append(_to_item(e))
        elif el.tag == "Models":
            doc_sum["Models"] = []
            for e in el:
                doc_sum["Models"].append(_to_item(e))
        elif el.tag == "Attributes":
            doc_sum["Attributes"] = []
            for e in el:
                doc_sum["Attributes"].append(_to_item(e))
        elif el.tag == "Status":
            doc_sum["Status"] = _to_item(el)
        elif el.tag == "Package":
            doc_sum["Package"] = _to_item(el)
            #for e in el:
            #    doc_sum["Package"].append(_to_item(e))
        elif el.tag == "Links":
            doc_sum["Links"] = []
            for _el in el:
                doc_sum["Links"].append(_to_item(_el))
            #    doc_sum["Package"].append(_to_item(e))
        elif el.tag == "Item":
            if "Name" in el.attrib and "Type" in el.attrib:
                if el.attrib["Type"] == "String":
                    doc_sum[el.attrib["Name"]] = str(el.text)
                elif el.attrib["Type"] == "Integer":
                    doc_sum[el.attrib["Name"]] = int(el.text)
                else:
                    print("attrb?", el, el.attrib, el.text)
            else:
                print("element?", el, el.attrib, el.text)
        else:
            print("tag?", el, el.attrib, el.text)
    return doc_sum


def parse_document_summary(tree):
    parse_function = {
        "RsUid": int,
        "GbUid": int,
        "AssemblyAccession": str,
        "LastMajorReleaseAccession": str,
        "LatestAccession": str,
        "ChainId": int,
        "AssemblyName": str,
        "UCSCName": str,
        "EnsemblName": str,
        "Taxid": int,
        "Organism": str,
        "SpeciesTaxid": int,
        "SpeciesName": str,
        "AssemblyType": str,
        "AssemblyStatus": str,
        "AssemblyStatusSort": str,
        "WGS": str,
        "Coverage": float,
        "ReleaseLevel": str,
        "ReleaseType": str,
        "SubmitterOrganization": str,
        "RefSeq_category": str,
        "FromType": str,
        "ContigN50": str,
        "ScaffoldN50": str,
        "FtpPath_GenBank": str,
        "FtpPath_RefSeq": str,
        "FtpPath_Assembly_rpt": str,
        "FtpPath_Stats_rpt": str,
        "FtpPath_Regions_rpt": str,
        "BioSampleAccn": str,
        "BioSampleId": str,
        "AsmReleaseDate_GenBank": str,
        "AsmReleaseDate_RefSeq": str,
        "SeqReleaseDate": str,
        "AsmUpdateDate": str,
        "SubmissionDate": str,
        "LastUpdateDate": str,
        "Primary": str,

        # bioproject
        "TaxId": int,
        "Project_Id": int,
        "Project_Acc": str,
        "Project_Type": str,
        "Project_Data_Type": str,
        "Sort_By_ProjectType": str,
        "Sort_By_DataType": str,
        "Sort_By_Organism": str,
        "Project_Target_Scope": str,
        "Project_Target_Material": str,
        "Project_Target_Capture": str,
        "Project_MethodType": str,
        "Registration_Date": str,
        "Project_Name": str,
        "Project_Title": str,
        "Organism_Name": str,
        "Organism_Strain": str,
        "Sequencing_Status": str,
        "Submitter_Organization": str,
        "Supergroup": str,
        "Relevance_Agricultural": str,
        "Relevance_Medical": str,
        "Relevance_Industrial": str,
        "Relevance_Environmental": str,
        "Relevance_Evolution": str,
        "Relevance_Model": str,
        "Relevance_Other": str,
    }
    doc = {}
    for el in tree:
        if el.tag in parse_function:
            try:
                doc[el.tag] = parse_function[el.tag](el.text)
            except Exception as ex:
                logger.warning(f'error: {el.tag} {el.text} {ex}')
        elif el.tag == "PropertyList":
            property_list = []
            for e in el:
                property_list.append(e.text)
            doc["PropertyList"] = property_list
        elif el.tag == "Synonym":
            synonym = {}
            for e in el:
                synonym[e.tag] = e.text
            doc["Synonym"] = synonym
        elif el.tag == "Meta":
            logger.debug(f'Meta: {el} {el.tag} {el.text} {el.attrib}')
        else:
            logger.debug(f'{el} {el.tag} {el.text} {el.attrib}')
    return doc


def parse_doc_sum(t):
    doc_sum = {}
    for el in t:
        if el.tag == "Id":
            doc_sum["id"] = el.text
        elif el.tag == "Item":
            if "Name" in el.attrib and "Type" in el.attrib:
                if el.attrib["Type"] == "String":
                    doc_sum[el.attrib["Name"]] = str(el.text)
                elif el.attrib["Type"] == "Integer":
                    doc_sum[el.attrib["Name"]] = int(el.text)
                else:
                    print("attrb?", el, el.attrib, el.text)
            else:
                print("element?", el, el.attrib, el.text)
        else:
            print("tag?", el, el.attrib, el.text)
    return doc_sum


class NcbiEutils:

    def __init__(self, api_key: str = None, cache_dir: str = None):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.api_key = api_key
        self.cache_dir = cache_dir

    @staticmethod
    def build_url_params(p):
        return "&".join(map(lambda x: str(x[0]) + "=" + str(x[1]), p.items()))

    def esearch(self, db, term, retmode="xml", retmax=10, retstart=0):
        params = {
            "retmode": retmode,
            "db": db,
            "term": term,
            "retmax": retmax,
            "retstart": retstart,
        }
        if self.api_key:
            params["api_key"] = self.api_key

        res = requests.get(
            self.base_url + "/esearch.fcgi?" + self.build_url_params(params)
        )

        tree = ET.fromstring(res.content)

        id_list = set()
        ret_max = None
        ret_start = None
        count = None
        query_translation = None
        for o in tree:
            if o.tag == "IdList":
                id_list = {i.text for i in o}
            elif o.tag == "Count":
                count = int(o.text)
            elif o.tag == "RetStart":
                ret_start = int(o.text)
            elif o.tag == "RetMax":
                ret_max = int(o.text)
            elif o.tag == "QueryTranslation":
                query_translation = o.text
            else:
                print(o.tag, o.text)

        return {
            "ids": id_list,
            "count": count,
            "ret_max": ret_max,
            "ret_start": ret_start,
            "query_translation": query_translation,
        }

    def _get_eutils_summary_content(self, db, record_id, ret_mode="xml"):
        if self.cache_dir:
            cache_file = f'{self.cache_dir}/{db}/{record_id}.xml'
            if os.path.exists(cache_file):
                with open(cache_file, 'rb') as fh:
                    return fh.read()

        params = {"retmode": ret_mode, "db": db, "id": record_id}
        if self.api_key:
            params["api_key"] = self.api_key

        res = requests.get(
            self.base_url + "/esummary.fcgi?" + self.build_url_params(params)
        )

        if self.cache_dir:
            cache_file = f'{self.cache_dir}/{db}/{record_id}.xml'
            with open(cache_file, 'w', encoding='utf-8') as fh:
                fh.write(res.content.decode('utf-8'))

        return res.content

    def _e_fetch(self, db, record_id, ret_mode="xml"):
        params = {"retmode": ret_mode, "db": db, "id": record_id}

        if self.api_key:
            params["api_key"] = self.api_key
        res = requests.get(
            self.base_url + "/efetch.fcgi?" + self.build_url_params(params)
        )

        return res.content

    def efetch(self):
        pass

    def esummary(self, db, record_id, ret_mode="xml"):
        content = self._get_eutils_summary_content(db, record_id, ret_mode)
        tree = ET.fromstring(content)

        docs = []
        for el in tree:
            if el.tag == "DocSum":
                docs.append(parse_doc_sum(el))
            elif el.tag == "DocumentSummarySet":
                docs.append(parse_document_summary_set(el))
            else:
                print(el.tag)
        return docs

    def elink(self, dbfrom, db, ids):
        params = {"dbfrom": dbfrom, "db": db, "id": ids}
        if self.api_key:
            params["api_key"] = self.api_key
        res = requests.get(
            self.base_url + "/elink.fcgi?" + self.build_url_params(params)
        )
        return res


class NcbiEutilsCacheMongo(NcbiEutils):

    DB_COLLECTION = {
        'assembly': 'ncbi_assembly',
        'bioproject': 'ncbi_bioproject',
        'nuccore': 'ncbi_nuccore',
        'protein': 'ncbi_protein',
    }

    def __init__(self, api_key, mongo_database):
        self.database = mongo_database
        super().__init__(api_key=api_key)

    def esummary(self, db, record_id, ret_mode="xml"):
        col_ncbi_assembly = self.database[NcbiEutilsCacheMongo.DB_COLLECTION[db]]

        logger.debug('search cache')

        doc = col_ncbi_assembly.find_one({'_id': record_id})

        if doc is None:
            logger.debug('fetch live')

            ret = super().esummary(db, record_id, ret_mode)
            ret = ret[0] if len(ret) > 1 else ret
            ret = ret[0] if len(ret) > 1 else ret

            print(ret)

            if len(ret) != 0:
                logger.debug('insert result in cache')
                ret['_id'] = record_id
                col_ncbi_assembly.insert_one(ret)

                doc = ret
            else:
                raise ValueError(f'invalid query {db}:{record_id}')

        return doc
