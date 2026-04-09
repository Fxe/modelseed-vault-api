import logging
import lxml.etree as ET


LOGGER = logging.getLogger(__name__)


"""
from modelseedpy_ext.re.etl.etl_transform_graph import ETLTransformGraph
old stuff
class TransformSbml(ETLTransformGraph):

    def __init__(self):
        super().__init__()

    def transform(self, sbml_file: str):
        nodes = {}
        edges = {}



transform = TransformSbml()
"""


def parse_sbml(elem, parser):
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'model':
            parse_model(elem, parser)
        elif action == 'end' and ET.QName(elem).localname == 'sbml':
            return None


def parse_model(elem, parser):
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'listOfUnitDefinitions':
            parse_list_of_unit_definitions(elem, parser)
        elif action == 'start' and ET.QName(elem).localname == 'listOfParameters':
            list_of_parameters = parse_list_of_parameters(elem, parser)
        elif action == 'start' and ET.QName(elem).localname == 'listOfCompartments':
            list_of_compartments = parse_list_of_compartments(elem, parser)
        elif action == 'start' and ET.QName(elem).localname == 'listOfSpecies':
            parse_list_of_species(elem, parser)
        elif action == 'start' and ET.QName(elem).localname == 'listOfReactions':
            parse_list_of_reactions(elem, parser)

        elif action == 'start' and ET.QName(elem).localname == 'listOfGeneProducts':
            parse_list_of_gene_products(elem, parser)
        elif action == 'start' and ET.QName(elem).localname == 'listOfObjectives':
            parse_list_of_objectives(elem, parser)

        elif action == 'end' and ET.QName(elem).localname == 'model':
            print(list_of_compartments, list_of_parameters)
            return None
        else:
            print('!', elem)


def parse_list_of_unit_definitions(elem, parser):
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'unitDefinition':
            pass
        elif action == 'end' and ET.QName(elem).localname == 'listOfUnitDefinitions':
            return None


def parse_list_of_compartments(elem, parser):
    result = []
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'compartment':
            o = dict(elem.attrib)
        elif action == 'end' and ET.QName(elem).localname == 'compartment':
            result.append(o)
        elif action == 'end' and ET.QName(elem).localname == 'listOfCompartments':
            return result
        else:
            LOGGER.warning(f'Ignore [listOfCompartments] tag: [{action}] {elem.tag}')


def parse_list_of_parameters(elem, parser):
    result = []
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'parameter':
            o = dict(elem.attrib)
        elif action == 'end' and ET.QName(elem).localname == 'parameter':
            result.append(o)
        elif action == 'end' and ET.QName(elem).localname == 'listOfParameters':
            return result
        else:
            LOGGER.warning(f'Ignore [listOfParameters] tag: [{action}] {elem.tag}')


def parse_species(elem, parser):
    pass


def parse_list_of_species(elem, parser):
    result = []
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'species':
            o = parse_species(elem, parser):
            elif action == 'end' and ET.QName(elem).localname == 'listOfSpecies':
            return result
        else:
            LOGGER.warning(f'Ignore [listOfSpecies] tag: [{action}] {elem.tag}')


def parse_list_of_reactions(elem, parser):
    result = []
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'reaction':
            pass
        elif action == 'end' and ET.QName(elem).localname == 'listOfReactions':
            return result
        else:
            LOGGER.warning(f'Ignore [listOfReactions] tag: [{action}] {elem.tag}')


def parse_list_of_objectives(elem, parser):
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'reaction':
            pass
        elif action == 'end' and ET.QName(elem).localname == 'listOfObjectives':
            return None


def parse_list_of_gene_products(elem, parser):
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'geneProduct':
            pass
        elif action == 'end' and ET.QName(elem).localname == 'listOfGeneProducts':
            return None


with open('./e_coli_core.xml', 'rb') as fh:
    parser = ET.iterparse(fh, events=('end', 'start'))
    for action, elem in parser:
        if action == 'start' and ET.QName(elem).localname == 'sbml':
            # print(elem)
            parse_sbml(elem, parser)
        else:
            print('error ')


class BasicElementParser:

    def __init__(self, chain, end):
        self.end = end
        self.chain = chain

    def parse(self, elem, parser):
        result = []
        for action, elem in parser:
            tag = ET.QName(elem).localname
            if action == 'start' and tag in self.chain:
                o = chain[tag](action, elem)
            if action == 'end' and tag in self.chain:
                result.append(o)
            elif action == 'end' and tag == self.end:
                return result
            else:
                LOGGER.warning(f'Ignore [{self.end}] tag: [{action}] {elem.tag}')

        raise Exception('parser failed to exit tag: ' + self.end)

## Load Model and Translate

def load_model(model_id, model_data, mapping, index=0, skip_gpr=False):
    from cobra.core import Model, Metabolite, Reaction
    from modelseedpy import MSBuilder

    metabolites_delete = set(mapping.get('metabolites_deleted', set()))

    metabolites = {}
    reactions = {}
    metabolite_translation = {}

    for m in model_data['metabolites']:
        if m['id'] not in metabolites_delete:
            m_comp = m['compartment']
            if m['compartment'] in mapping['compartments']:
                m_comp = mapping['compartments'][m['compartment']]
            m_id = m['id']
            if m['id'] in mapping['metabolites']:
                t = mapping['metabolites'][m['id']]
                m_id = f'{t}_{m_comp}{index}'
                metabolite_translation[m['id']] = m_id
            metabolite = Metabolite(m_id, m.get('formula'), m['name'], 0, m_comp)
            if metabolite.id in metabolites:
                print('duplicate')
            metabolites[metabolite.id] = metabolite

    for r in model_data['reactions']:
        if r['id'] not in mapping['exchanges']:
            reaction = Reaction(r['id'], r['name'], r.get('subsystem'), r['lower_bound'], r['upper_bound'])
            if not skip_gpr:
                reaction.gene_reaction_rule = r['gene_reaction_rule']
            reaction_metabolites = {metabolites[metabolite_translation.get(m, m)]: v for m, v in
                                    r['metabolites'].items()}
            reaction.add_metabolites(reaction_metabolites)
            if reaction.id in reactions:
                print('duplicate')
            reactions[reaction.id] = reaction

    model = Model(model_id)
    model.add_metabolites(list(metabolites.values()))
    model.add_reactions(list(reactions.values()))
    MSBuilder.add_exchanges_to_model(model, extra_cell=f'e')

    return model


def example():
    from cobra.core import Model, Metabolite, Reaction
    import json
    model_load_data = {
        'iAF692': ['/home/fliu/scratch/data/sbml/iAF692/iAF692.json',
                   '/home/fliu/scratch/data/sbml/iAF692/mapping.json'],
        'iMB745': ['/home/fliu/scratch/data/sbml/iMB745/iMB745.json',
                   '/home/fliu/scratch/data/sbml/iMB745/mapping.json'],
        'iVS941': ['/home/fliu/scratch/data/sbml/iVS941/iVS941.json',
                   '/home/fliu/scratch/data/sbml/iVS941/mapping.json'],
        'iMAC868': ['/home/fliu/scratch/data/sbml/iMAC868/iMAC868.json',
                    '/home/fliu/scratch/data/sbml/iMAC868/mapping.json'],
        # 'NmrFL413': ['/home/fliu/scratch/data/sbml/NmrFL413/NmrFL413.json', '/home/fliu/scratch/data/sbml/NmrFL413/mapping.json'],
        'iMN22HE': ['/home/fliu/scratch/data/sbml/iMN22HE/iMN22HE.json',
                    '/home/fliu/scratch/data/sbml/iMN22HE/mapping.json'],
    }

    models = {}

    for model_id, (file_model, file_mapping) in model_load_data.items():
        print(model_id)
        model_data = None
        mapping = None
        with open(file_model, 'r') as fh:
            model_data = json.load(fh)
        with open(file_mapping, 'r') as fh:
            mapping = json.load(fh)

        model_translated = load_model(model_id, model_data, mapping, index='', skip_gpr=True)
        obj = Reaction('obj_atp', 'Objective ATP', '', 0, 1000)
        obj.add_metabolites({
            model_translated.metabolites.cpd00067_c: 1,
            model_translated.metabolites.cpd00008_c: 1,
            model_translated.metabolites.cpd00009_c: 1,
            model_translated.metabolites.cpd00002_c: -1,
            model_translated.metabolites.cpd00001_c: -1
        })
        model_translated.add_reactions([obj])
        model_translated.objective = 'obj_atp'
        models[model_id] = model_translated

    import pandas as pd
    from modelseedpy.core.msmedia import MSMedia
    from modelseedpy_ext.utils import load_medias
    def load_medias(filep, index_col=0):
        medias = pd.read_csv(filep, sep='\t', index_col=index_col).to_dict()
        return medias

    medias = load_medias('/home/fliu/python3/ModelSEEDpy-Ext/modelseedpy_ext/profiler/medias_core.tsv', 2)
    media_ids = [
        'Glc/O2', 'Ac/O2', 'Fum/O2', 'Succ/O2', 'Akg/O2', 'Cit/O2', 'LMal', 'LLac/O2', 'Dlac/O2',
        'Pyr/O2', 'Glyc/O2', 'Etho/O2', 'For/O2',
        'Glc', 'Ac', 'Fum', 'Succ', 'Akg', 'Cit', 'mal-L', 'Llac', 'Dlac',
        'Pyr', 'Glyc', 'Etho', 'For',
        'For/NO2', 'For/NO3', 'For/NO',
        'Pyr/NO2', 'Pyr/NO3', 'Pyr/NO',
        'Ac/NO2', 'Ac/NO3', 'Ac/NO',
        'Glc/DMSO', 'Glc/TMAO', 'Pyr/DMSO', 'Pyr/TMAO',
        'Pyr/SO4', 'Pyr/SO3',
        'H2/CO2', 'H2/Ac', 'For/SO4/H2', 'LLac/SO4/H2', 'For/SO4', 'LLac/SO4',
        'H2/SO4',
        'empty', 'Light', 'ANME', 'Methane'
    ]
    media_data = {}
    for media_id in media_ids:
        d = medias[media_id]
        media_const = {
            'cpd00001': 1000,
            'cpd00067': 1000
        }
        for cpd_id in d:
            v = d[cpd_id]
            if v != 0:
                media_const[cpd_id] = v
        media = MSMedia.from_dict(media_const)
        media_data[media_id] = media

    class AtpCoreProfiler:

        def __init__(self, master, medias, atp_hydrolysis_id, extracell='e0', complement=None):
            self.master = master
            self.atp_hydrolysis_id = atp_hydrolysis_id
            self.medias = medias
            if complement is None:
                self.comp = {'EX_cpd00001_e0': 1000, 'EX_cpd00067_e0': 1000}
            else:
                self.comp = complement
            self.extracell = extracell

        def profile_genome(self, genome_id):
            model = self.master.fn_get_model(genome_id)
            model.objective = self.atp_hydrolysis_id
            return self.test_core_model(model, self.medias, cmp=self.extracell)

        def test_media(self, model, media, exs, cmp):
            media_const = media.get_media_constraints(cmp)
            medium_ = {}
            for cpd_id in media_const:
                if cpd_id not in exs:
                    return None
            for cpd_id in media_const:
                lb, ub = media_const[cpd_id]
                if cpd_id in exs:
                    rxn_exchange = exs[cpd_id]
                    medium_[rxn_exchange.id] = -1 * lb
            medium_.update(self.comp)
            model.medium = medium_
            # print(media.id, model.medium)
            sol = model.optimize()
            return sol

        def get_metabolite_exchanges(self, model):
            metabolite_exchanges = {}
            for rxn_exchange in model.exchanges:
                metabolites = rxn_exchange.metabolites
                if len(metabolites) == 1:
                    metabolite_exchanges[list(metabolites)[0].id] = rxn_exchange
                else:
                    print('ignore', rxn_exchange)
            return metabolite_exchanges

        def test_core_model(self, model, medias, v=False, cmp='e0'):
            res = {}
            media_out = {}
            metabolite_exchanges = self.get_metabolite_exchanges(model)
            for media_id in medias:
                # print(media_id)
                sol = self.test_media(model, medias[media_id], metabolite_exchanges, cmp)
                media_out[media_id] = None
                res[media_id] = 0
                if sol and sol.status == 'optimal':
                    ex_out = dict(filter(lambda x: x[0].startswith('EX_') and x[1] != 0, sol.fluxes.to_dict().items()))
                    media_out[media_id] = ex_out
                    # print(media_id, sol.objective_value, ex_out)
                    atp_val = sol.objective_value
                    if atp_val > 0:
                        for cpd_id in medias[media_id].get_media_constraints(cmp):
                            media_cpd_val = sol.fluxes[metabolite_exchanges[cpd_id].id]
                            if media_cpd_val == 0:
                                atp_val = 0
                                # print('reject', media_id)
                    res[media_id] = atp_val
            return res, media_out

    class GeProf:

        def __init__(self, genome_ids, fn_get_model, fn_get_genome, genome_order=None):
            self.genome_order = None
            self.genome_ids = genome_ids
            self.fn_get_model = fn_get_model
            self.fn_get_genome = fn_get_genome

        def profile(self):
            pass

    ge_prof = GeProf(['iAF692', 'iMB745', 'iVS941', 'iMAC868', 'iMN22HE'], lambda x: models[x], None)
    profiler = AtpCoreProfiler(ge_prof, media_data, 'obj_atp', 'e', {})
    for model in models.values():
        if 'ATPM' in model.reactions:
            r = model.reactions.get_by_id('ATPM')
            r.lower_bound = 0
            r.upper_bound = 0

    profiler.profile_genome('iMAC868')

    medium_iMB745 = {
        'EX_cpd00001_e': 1000,
        'EX_4hphac_e': 1,
        'EX_cpd00029_e': 1,
        'EX_actn_R_e': 1,
        'EX_cpd00035_e': 1,
        'EX_alac_S_e': 1,
        'EX_cpd01024_e': 1,
        'EX_cpd00011_e': 1,
        'EX_cpd00204_e': 1,
        'EX_cobalt2_e': 1,
        'EX_cpd00084_e': 1,
        'EX_cpd00425_e': 1,
        'EX_dms_e': 1,
        'EX_cpd10515_e': 1,
        'EX_cpd10516_e': 1,
        'EX_cpd00047_e': 1,
        'EX_cpd00023_e': 1,
        'EX_cpd00033_e': 1,
        'EX_cpd11640_e': 1,

        'EX_cpd00239_e': 1,
        'EX_cpd00067_e': 1000,
        'EX_hco3_e': 1,
        'EX_cpd00322_e': 1,
        'EX_cpd00107_e': 1,
        'EX_cpd00039_e': 1,
        'EX_meoh_e': 1,
        'EX_mg2_e': 1,
        'EX_cpd00187_e': 1,
        'EX_mn2_e': 1,
        'EX_mobd_e': 1,
        'EX_cpd00528_e': 1,
        'EX_cpd00971_e': 1,
        'EX_nac_e': 1,
        'EX_cpd00013_e': 1,
        'EX_cpd00009_e': 1,
        'EX_cpd00644_e': 1,
        'EX_cpd00129_e': 1,
        'EX_cpd00020_e': 1,
        'EX_ribflv_e': 1,
        'EX_cpd00048_e': 1,
        'EX_cpd00305_e': 1,
        'EX_cpd00441_e': 1,
        'EX_trp_L_e': 1,
        'EX_unknown_cbl1deg_e': 1,
        'EX_unknown_rbfdeg_e': 1,
        'EX_cpd00156_e': 1,
        'EX_wo4_e': 1,
        'EX_zn2_e': 1}

    medium_iVS941 = {
        'EX_cpd00001_e': 1000,
        'EX_cpd00067_e': 1000,
        'EX_4abz_e': 1,
        'EX_4hphac_e': 1,
        'EX_CO2_e': 1000,
        'EX_H2CO3_e': 1000,
        'EX_cpd00029_e': 1,
        'EX_actn_R_e': 1,
        'EX_cpd00035_e': 1,
        'EX_alac_S_e': 1,
        'EX_btn_e': 1000,
        'EX_ca2_e': 1000,
        'EX_cbi_e': 1000,
        'EX_cbl1_e': 1000,
        'EX_cbl1hbi_e': 1000,
        'EX_cd2_e': 1000,
        'EX_cpd01024_e': 1000,
        'EX_ch4s_e': 1000,
        'EX_cit_e': 1000,
        'EX_cl_e': 1000,
        'EX_cpd00011_e': 1000,
        'EX_cpd00204_e': 1000,
        'EX_cobalt2_e': 1000,
        'EX_cu2_e': 1000,
        'EX_cpd00084_e': 1000,
        'EX_cpd00425_e': 1000,
        'EX_dms_e': 1000,
        'EX_etha_e': 1000,
        'EX_cpd10515_e': 1000,
        'EX_cpd10516_e': 1000,
        'EX_fol_e': 1000,
        'EX_gcald_e': 1000,
        'EX_glcn_e': 1000,
        'EX_cpd00023_e': 1000,
        'EX_cpd00033_e': 1000,
        'EX_glyb_e': 1000,
        'EX_cpd11640_e': 1000,

        'EX_cpd00239_e': 1000,

        'EX_cpd00322_e': 1000,
        'EX_ind3ac_e': 1000,
        'EX_k_e': 1000,
        'EX_cpd00107_e': 1000,
        'EX_cpd00039_e': 1000,
        'EX_meoh_e': 1000,
        'EX_mg2_e': 1000,
        'EX_cpd00187_e': 1000,
        'EX_mn2_e': 1000,
        'EX_mphen_e': 1000,
        'EX_mphenh2_e': 1000,
        'EX_cpd00528_e': 1000,
        'EX_cpd00971_e': 1000,
        'EX_nac_e': 1000,
        'EX_cpd00013_e': 1000,
        'EX_ni2_e': 1000,
        'EX_pac_e': 1000,
        'EX_cpd00009_e': 1000,
        'EX_cpd00644_e': 1000,
        'EX_cpd00129_e': 1000,
        'EX_cpd00020_e': 1000,
        'EX_ribflv_e': 1000,
        'EX_s_e': 1000,
        'EX_so3_e': 1000,
        'EX_cpd00048_e': 1000,
        'EX_cpd00305_e': 1000,
        'EX_cpd00441_e': 1000,
        'EX_tsul_e': 1000,
        'EX_unknown_cbl1deg_e': 1000,
        'EX_unknown_rbfdeg_e': 1000,
        'EX_cpd00156_e': 1000,
        'EX_zn2_e': 1000}
    model.objective = 'obj_atp'
    model.medium = medium_iVS941
    model.summary()

    solution = model.optimize()
    for r in model.reactions:
        v = solution[r.id]
        # if v != 0:
        #    print(v, r)
        if v > 1e-6 or v < -1e-6:
            print(v, r)

    model.reactions.MRP.lower_bound = 0
    model.reactions.MRP.upper_bound = 0
    model.metabolites.cpd00067_c.summary(solution)

    model.objective = model.reactions.overall.id

    model.medium = {'EX_4abz_e': 1000,
                    'EX_4hphac_e': 1000,
                    'EX_cpd00029_e': 1000,
                    'EX_actn_R_e': 1000,
                    'EX_cpd00035_e': 1000,
                    'EX_alac_S_e': 1000,
                    'EX_biomass_met_e': 1000,
                    'EX_btn_e': 1000,
                    'EX_ca2_e': 1000,
                    'EX_cbi_e': 1000,
                    'EX_cbl1_e': 1000,
                    'EX_cbl1hbi_e': 1000,
                    'EX_cd2_e': 1000,
                    'EX_cpd01024_e': 1000,
                    'EX_ch4s_e': 1000,
                    'EX_cit_e': 1000,
                    'EX_cl_e': 1000,
                    'EX_cpd00011_e': 1000,
                    'EX_cpd00204_e': 1000,
                    'EX_cobalt2_e': 1000,
                    'EX_cu2_e': 1000,
                    'EX_cpd00084_e': 0,
                    'EX_cpd00425_e': 0,
                    'EX_dms_e': 1000,
                    'EX_etha_e': 1000,
                    'EX_cpd10515_e': 1000,
                    'EX_cpd10516_e': 1000,
                    'EX_fol_e': 1000,
                    'EX_cpd00047_e': 1000,
                    'EX_gcald_e': 1000,
                    'EX_glcn_e': 1000,
                    'EX_cpd00023_e': 1000,
                    'EX_cpd00033_e': 1000,
                    'EX_glyb_e': 1000,
                    'EX_cpd11640_e': 1000,
                    'EX_cpd00001_e': 1000,
                    'EX_cpd00239_e': 1000,
                    'EX_cpd00067_e': 1000,
                    'EX_hco3_e': 1000,
                    'EX_cpd00322_e': 1000,
                    'EX_ind3ac_e': 1000,
                    'EX_k_e': 1000,
                    'EX_cpd00107_e': 1000,
                    'EX_cpd00039_e': 1000,
                    'EX_meoh_e': 0,
                    'EX_mg2_e': 1000,
                    'EX_cpd00187_e': 1000,
                    'EX_mn2_e': 1000,
                    'EX_mobd_e': 1000,
                    'EX_cpd00528_e': 1000,
                    'EX_cpd00971_e': 1000,
                    'EX_nac_e': 1000,
                    'EX_cpd00013_e': 1000,
                    'EX_ni2_e': 1000,
                    'EX_pac_e': 1000,
                    'EX_cpd00009_e': 1000,
                    'EX_cpd00644_e': 1000,
                    'EX_cpd00129_e': 1000,
                    'EX_cpd00020_e': 0,
                    'EX_ribflv_e': 1000,
                    'EX_cpd00048_e': 1000,
                    'EX_cpd00305_e': 1000,
                    'EX_cpd00441_e': 10,
                    'EX_trp_L_e': 1000,
                    'EX_unknown_cbl1deg_e': 1000,
                    'EX_unknown_rbfdeg_e': 1000,
                    'EX_cpd00156_e': 1000,
                    'EX_wo4_e': 1000,
                    'EX_zn2_e': 1000}
    model.summary()
    for r in model.reactions:
        if r.lower_bound != 0 and r.lower_bound != -1000:
            print(r.id, r.lower_bound, r.upper_bound)
        if r.upper_bound != 0 and r.upper_bound != 1000:
            print(r.id, r.lower_bound, r.upper_bound)
    profiler.profile_genome('iMB745')

    profiler.profile_genome('iAF692')