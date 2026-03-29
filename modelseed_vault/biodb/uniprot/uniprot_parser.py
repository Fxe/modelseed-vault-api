import logging
import lxml.etree as ET

logging.getLogger(__name__)


class Wut:
    def __init__(self, tag, verbose=False):
        self.tag = tag
        self.ns = ""
        self.ns_tag = tag
        self.verbose = verbose

    def parse_section(self, elem, parser, tag_end, capture):
        # print('section', tag_end, capture)
        res = dict(elem.attrib)
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.ns + tag_end:
                return res
            elif action in capture and tag in capture[action]:
                fn, t, e = capture[action][tag]
                if t == list:
                    if e not in res:
                        res[e] = []
                    res[e].append(fn(elem, parser))
            else:
                if self.verbose:
                    print(
                        f"[no capture] [section:{tag_end}] [tag:{self.tag}]",
                        action,
                        tag,
                        elem.attrib,
                        elem.text,
                        elem.sourceline,
                    )

    def parse_blob(self, elem, parser, tag_end):
        res = dict(elem.attrib)
        if elem.text and len(elem.text.strip()) > 0:
            res["value"] = elem.text
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.ns + tag_end:
                # sometimes text is captured at the end tag !?
                if elem.text:
                    res["value"] = elem.text
                return res
            else:
                if self.verbose:
                    print(
                        f"[no capture parse_blob] {self.tag}",
                        action,
                        tag,
                        elem.attrib,
                        elem.text,
                        elem.sourceline,
                    )

    def parse(self, elem, parser):
        # print('!!!!!!!')
        res = dict(elem.attrib)
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.ns_tag:
                return res
            else:
                if self.verbose:
                    print(
                        f"[no capture] {self.tag}", action, tag, elem.attrib, elem.text
                    )


class CaptureDbReference(Wut):
    def parse(self, elem, parser):
        res = dict(elem.attrib)
        res["xml_sourceline"] = elem.sourceline
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.ns_tag:
                return res
            elif action == "start" and tag == self.ns + "property":
                if "property" not in res:
                    res["property"] = []
                res["property"].append(self.parse_blob(elem, parser, "property"))
            elif action == "start" and tag == self.ns + "molecule":
                if "molecule" not in res:
                    res["molecule"] = []
                res["molecule"].append(self.parse_blob(elem, parser, "molecule"))
            else:
                if self.verbose:
                    print(
                        f"[no capture] [tag:{self.tag}]",
                        action,
                        tag,
                        elem.attrib,
                        elem.text,
                        elem.sourceline,
                    )


class CaptureEvidence(Wut):
    def parse_source(self, elem, parser):
        res = dict(elem.attrib)
        res["dbReference"] = []
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.ns + "source":
                return res
            elif action == "start" and tag == self.ns + "dbReference":
                blob = dict(elem.attrib)
                blob["value"] = elem.text
                res["dbReference"].append(blob)
            elif action == "end" and tag == self.ns + "dbReference":
                pass
            else:
                print(
                    f"[no capture source] {self.tag}",
                    action,
                    tag,
                    elem.attrib,
                    elem.text,
                )

    def parse(self, elem, parser):
        res = dict(elem.attrib)
        res["source"] = []
        res["importedFrom"] = []
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and tag == self.ns + "source":
                res["source"].append(self.parse_source(elem, parser))
            elif action == "end" and tag == self.ns_tag:
                return res
            else:
                if self.verbose:
                    print(
                        f"[no capture] {self.tag}", action, tag, elem.attrib, elem.text
                    )


class CaptureDbReferenceType(Wut):

    pass


class CaptureProtein(Wut):
    def parse_name(self, elem, parser, end):
        res = dict(elem.attrib)
        for _action, _elem in parser:
            _tag = _elem.tag
            if _action == "end" and _tag == self.ns + end:
                return res
            elif _action == "start" and _tag == self.ns + "fullName":
                o = dict(_elem.attrib)
                o["value"] = _elem.text
                res["full"] = o
            elif _action == "end" and _tag == self.ns + "fullName":
                pass
            elif _action == "start" and _tag == self.ns + "shortName":
                o = dict(_elem.attrib)
                o["value"] = _elem.text
                res["short"] = o
            elif _action == "end" and _tag == self.ns + "shortName":
                pass
            elif _action == "start" and _tag == self.ns + "ecNumber":
                o = dict(_elem.attrib)
                o["value"] = _elem.text
                res["ec"] = o
            elif _action == "end" and _tag == self.ns + "ecNumber":
                pass
            else:
                print(
                    f"[no _capture name] {_tag}",
                    _action,
                    _tag,
                    _elem.attrib,
                    _elem.text,
                )

    def parse_component(self, elem, parser):
        res = dict(elem.attrib)
        res["recommended_name"] = []
        res["alternative_name"] = []
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.ns + "component":
                return res
            elif action == "start" and tag == self.ns + "recommendedName":
                res["recommended_name"].append(
                    self.parse_name(elem, parser, "recommendedName")
                )
            elif action == "start" and tag == self.ns + "alternativeName":
                res["alternative_name"].append(
                    self.parse_name(elem, parser, "alternativeName")
                )
            else:
                if self.verbose:
                    print(
                        f"[no capture component] {self.tag}",
                        action,
                        tag,
                        elem.attrib,
                        elem.text,
                    )

    def parse(self, elem, parser):
        # print('!!!!!!!')
        res = dict(elem.attrib)
        res["recommended_name"] = []
        res["alternative_name"] = []
        res["component"] = []
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and tag == self.ns + "recommendedName":
                res["recommended_name"].append(
                    self.parse_name(elem, parser, "recommendedName")
                )
            elif action == "start" and tag == self.ns + "alternativeName":
                res["alternative_name"].append(
                    self.parse_name(elem, parser, "alternativeName")
                )
            elif action == "start" and tag == self.ns + "component":
                res["component"].append(self.parse_component(elem, parser))
            elif action == "end" and tag == self.ns_tag:
                return res
            else:
                if self.verbose:
                    print(
                        f"[no capture] {self.tag}", action, tag, elem.attrib, elem.text
                    )


class CaptureComment(Wut):
    def parse_isoform(self, elem, parser):
        res = dict(elem.attrib)
        res["id"] = []
        res["name"] = []
        res["sequence"] = []

        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.ns + "isoform":
                return res
            elif action == "start" and tag == self.ns + "id":
                blob = dict(elem.attrib)
                blob["value"] = elem.text
                res["id"].append(blob)
            elif action == "end" and tag == self.ns + "id":
                pass
            elif action == "start" and tag == self.ns + "name":
                blob = dict(elem.attrib)
                blob["value"] = elem.text
                res["name"].append(blob)
            elif action == "end" and tag == self.ns + "name":
                pass
            elif action == "start" and tag == self.ns + "sequence":
                blob = dict(elem.attrib)
                blob["value"] = elem.text
                res["sequence"].append(blob)
            elif action == "end" and tag == self.ns + "sequence":
                pass
            else:
                if self.verbose:
                    print(
                        f"[no capture parse_isoform] {self.tag}",
                        action,
                        tag,
                        elem.attrib,
                        elem.text,
                        elem.sourceline,
                    )

    def parse_conflict(self, elem, parser):
        res = dict(elem.attrib)
        res = {"sequence": []}
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.ns + "conflict":
                return res
            elif action == "start" and tag == self.ns + "sequence":
                blob = dict(elem.attrib)
                blob["value"] = elem.text
                res["sequence"].append(blob)
            elif action == "end" and tag == self.ns + "sequence":
                pass
            else:
                if self.verbose:
                    print(
                        f"[no capture parse_conflict] {self.tag}",
                        action,
                        tag,
                        elem.attrib,
                        elem.text,
                        elem.sourceline,
                    )

    def parse_location_blob(self, elem, parser):
        return self.parse_blob(elem, parser, "location")

    def parse_topology(self, elem, parser):
        return self.parse_blob(elem, parser, "topology")

    def parse_db_reference(self, elem, parser):
        return self.parse_blob(elem, parser, "dbReference")

    def parse_acronym(self, elem, parser):
        return self.parse_blob(elem, parser, "acronym")

    def parse_name(self, elem, parser):
        return self.parse_blob(elem, parser, "name")

    def parse_id(self, elem, parser):
        return self.parse_blob(elem, parser, "id")

    def parse_event(self, elem, parser):
        return self.parse_blob(elem, parser, "event")

    def parse_experiments(self, elem, parser):
        return self.parse_blob(elem, parser, "experiments")

    def parse_label(self, elem, parser):
        return self.parse_blob(elem, parser, "label")

    def parse_molecule(self, elem, parser):
        return self.parse_blob(elem, parser, "molecule")

    def parse_description(self, elem, parser):
        return self.parse_blob(elem, parser, "description")

    def parse_text(self, elem, parser):
        return self.parse_blob(elem, parser, "text")

    def parse_km(self, elem, parser):
        return self.parse_blob(elem, parser, "KM")

    def parse_vmax(self, elem, parser):
        return self.parse_blob(elem, parser, "Vmax")

    def parse_max(self, elem, parser):
        return self.parse_blob(elem, parser, "max")

    def parse_orientation(self, elem, parser):
        return self.parse_blob(elem, parser, "orientation")

    def parse_position(self, elem, parser):
        return self.parse_blob(elem, parser, "position")

    def parse_link(self, elem, parser):
        return self.parse_blob(elem, parser, "link")

    def parse_organisms_differ(self, elem, parser):
        return self.parse_blob(elem, parser, "organismsDiffer")

    def parse_disease(self, elem, parser):
        capture = {
            "start": {
                self.ns + "name": (self.parse_name, list, "name"),
                self.ns + "acronym": (self.parse_acronym, list, "acronym"),
                self.ns + "description": (self.parse_description, list, "description"),
                self.ns + "dbReference": (self.parse_db_reference, list, "dbReference"),
            }
        }
        return self.parse_section(elem, parser, "disease", capture)

    def parse_interactant(self, elem, parser):
        capture = {
            "start": {
                self.ns + "id": (self.parse_id, list, "id"),
                self.ns + "label": (self.parse_label, list, "label"),
                self.ns + "dbReference": (self.parse_db_reference, list, "dbReference"),
            }
        }
        return self.parse_section(elem, parser, "interactant", capture)

    def parse_reaction(self, elem, parser):
        capture = {
            "start": {
                self.ns + "text": (self.parse_text, list, "text"),
                self.ns + "dbReference": (self.parse_db_reference, list, "dbReference"),
            }
        }
        return self.parse_section(elem, parser, "reaction", capture)

    def parse_physiological_reaction(self, elem, parser):
        capture = {
            "start": {
                self.ns + "dbReference": (self.parse_db_reference, list, "dbReference")
            }
        }
        return self.parse_section(elem, parser, "physiologicalReaction", capture)

    def parse_kinetics(self, elem, parser):
        capture = {
            "start": {
                self.ns + "KM": (self.parse_km, list, "KM"),
                self.ns + "text": (self.parse_text, list, "text"),
                self.ns + "Vmax": (self.parse_vmax, list, "Vmax"),
            }
        }
        return self.parse_section(elem, parser, "kinetics", capture)

    def parse_cofactor(self, elem, parser):
        capture = {
            "start": {
                self.ns + "name": (self.parse_name, list, "name"),
                self.ns + "dbReference": (self.parse_db_reference, list, "dbReference"),
            }
        }
        return self.parse_section(elem, parser, "cofactor", capture)

    def parse_subcellular_location(self, elem, parser):
        capture = {
            "start": {
                self.ns + "location": (self.parse_location_blob, list, "location"),
                self.ns + "topology": (self.parse_topology, list, "topology"),
                self.ns + "orientation": (self.parse_orientation, list, "orientation"),
            }
        }
        return self.parse_section(elem, parser, "subcellularLocation", capture)

    def parse_ph(self, elem, parser):
        capture = {"start": {self.ns + "text": (self.parse_text, list, "text")}}
        return self.parse_section(elem, parser, "phDependence", capture)

    def parse_temperature(self, elem, parser):
        capture = {"start": {self.ns + "text": (self.parse_text, list, "text")}}
        return self.parse_section(elem, parser, "temperatureDependence", capture)

    def parse_absorption(self, elem, parser):
        capture = {
            "start": {
                self.ns + "text": (self.parse_text, list, "text"),
                self.ns + "max": (self.parse_max, list, "max"),
            }
        }
        return self.parse_section(elem, parser, "absorption", capture)

    def parse_location_section(self, elem, parser):
        capture = {
            "start": {self.ns + "position": (self.parse_position, list, "position")}
        }
        return self.parse_section(elem, parser, "location", capture)

    def parse_redox(self, elem, parser):
        capture = {"start": {self.ns + "text": (self.parse_text, list, "text")}}
        return self.parse_section(elem, parser, "redoxPotential", capture)

    def parse(self, elem, parser):
        res = dict(elem.attrib)
        res["xml_sourceline"] = elem.sourceline
        capture = {
            "start": {
                self.ns + "phDependence": (self.parse_ph, list, "phDependence"),
                self.ns
                + "temperatureDependence": (
                    self.parse_temperature,
                    list,
                    "temperatureDependence",
                ),
                self.ns + "cofactor": (self.parse_cofactor, list, "cofactor"),
                self.ns + "kinetics": (self.parse_kinetics, list, "kinetics"),
                self.ns
                + "subcellularLocation": (
                    self.parse_subcellular_location,
                    list,
                    "subcellularLocation",
                ),
                self.ns + "reaction": (self.parse_reaction, list, "reaction"),
                self.ns
                + "physiologicalReaction": (
                    self.parse_physiological_reaction,
                    list,
                    "physiologicalReaction",
                ),
                self.ns + "disease": (self.parse_disease, list, "disease"),
                self.ns + "interactant": (self.parse_interactant, list, "interactant"),
                self.ns + "absorption": (self.parse_absorption, list, "absorption"),
                self.ns + "location": (self.parse_location_section, list, "location"),
                self.ns + "redoxPotential": (self.parse_redox, list, "redoxPotential"),
                self.ns + "conflict": (self.parse_conflict, list, "conflict"),
                self.ns + "isoform": (self.parse_isoform, list, "isoform"),
                self.ns + "link": (self.parse_link, list, "link"),
                self.ns + "text": (self.parse_text, list, "text"),
                self.ns + "event": (self.parse_event, list, "event"),
                self.ns + "molecule": (self.parse_molecule, list, "molecule"),
                self.ns + "experiments": (self.parse_experiments, list, "experiments"),
                self.ns
                + "organismsDiffer": (
                    self.parse_organisms_differ,
                    list,
                    "organismsDiffer",
                ),
            }
        }
        for action, elem in parser:
            tag = elem.tag
            if action in capture and tag in capture[action]:
                fn, t, e = capture[action][tag]
                if t == list:
                    if e not in res:
                        res[e] = []
                    res[e].append(fn(elem, parser))
            elif action == "end" and tag == self.ns_tag:
                return res
            else:
                if self.verbose:
                    print(
                        f"[no capture] {self.tag}",
                        action,
                        tag,
                        elem.attrib,
                        elem.text,
                        elem.sourceline,
                    )


class CaptureGene(Wut):
    def parse(self, elem, parser):
        # print('!!!!!!!')
        res = dict(elem.attrib)
        res["name"] = []
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and tag == self.ns + "name":
                blob = dict(elem.attrib)
                blob["value"] = elem.text
                res["name"].append(blob)
            elif action == "end" and tag == self.ns + "name":
                pass
            elif action == "end" and tag == self.ns_tag:
                return res
            else:
                if self.verbose:
                    print(
                        f"[no capture] {self.tag}", action, tag, elem.attrib, elem.text
                    )


class SwissProtParser:
    def __init__(self, parse_limit=None):
        """

        :param parse_limit: max number of entries
        """
        self.max_entries = parse_limit

        self.tag_entry = "entry"
        self.tag_accession = "accession"
        self.tag_name = "name"
        self.tag_sequence = "sequence"

        self.tag_capture = {
            "protein": (CaptureProtein("protein", False), None),
            "gene": (CaptureGene("gene", True), list),
            "organism": (Wut("organism"), None),
            "organismHost": (Wut("organismHost"), list),
            "geneLocation": (Wut("geneLocation"), list),
            "reference": (Wut("reference"), list),
            "comment": (CaptureComment("comment", True), list),
            "dbReference": (CaptureDbReference("dbReference", True), list),
            "proteinExistence": (Wut("proteinExistence", True), None),
            "keyword": (Wut("keyword"), list),
            "feature": (Wut("feature"), list),
            "evidence": (CaptureEvidence("evidence", True), list),
        }
        self.ns_map = {}
        # self.tag_protein_existence = 'proteinExistence'

    def set_xml_ns(self, ns):
        for x in self.tag_capture:
            self.ns_map[ns + x] = x
            self.tag_capture[x][0].ns_tag = ns + x
            self.tag_capture[x][0].ns = ns
        self.ns_map = dict({(ns + x, x) for x in self.tag_capture})
        self.tag_entry = ns + self.tag_entry
        self.tag_accession = ns + self.tag_accession
        self.tag_name = ns + self.tag_name
        self.tag_sequence = ns + self.tag_sequence

    def parse(self, fh):
        yielded = 0
        parser = ET.iterparse(fh, events=("end", "start"))
        for action, elem in parser:
            tag = elem.tag
            logging.debug(f'parse {action} {tag}')
            if action == "start" and tag == self.tag_entry:
                yield self.parse_entry(elem, parser)
                yielded += 1
            else:
                logging.debug(f'parse {action} {elem}')
            if self.max_entries and yielded >= self.max_entries:
                break

    def parse_reference(self, elem, parser, end):
        # print('!!!!!!!')
        res = dict(elem.attrib)
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == end:
                return res
            else:
                # print('[no capture] parse_reference', action, tag, elem.attrib, elem.text)
                pass

    def parse_entry(self, elem, parser):

        uniprot_entry = dict(elem.attrib)
        uniprot_entry["xml_sourceline"] = elem.sourceline
        uniprot_entry["name"] = []
        uniprot_entry["accession"] = []

        for x in self.tag_capture:
            if x not in uniprot_entry:
                if self.tag_capture[x][1] and self.tag_capture[x][1] == list:
                    uniprot_entry[x] = []

        prev_text = None
        for action, elem in parser:
            tag = elem.tag

            if action == "start" and tag in self.ns_map:
                tag_base = self.ns_map[tag]
                res = self.tag_capture[tag_base][0].parse(elem, parser)
                if (
                    self.tag_capture[tag_base][1]
                    and self.tag_capture[tag_base][1] == list
                ):
                    uniprot_entry[tag_base].append(res)
                else:
                    uniprot_entry[tag_base] = res

            elif action == "start" and tag == self.tag_sequence:
                uniprot_entry["protein_sequence"] = {"value": elem.text}
                uniprot_entry["protein_sequence"].update(elem.attrib)
            elif action == "end" and tag == self.tag_sequence:
                if uniprot_entry["protein_sequence"]["value"] is None:
                    uniprot_entry["protein_sequence"]["value"] = elem.text
            elif action == "end" and tag == self.tag_entry:
                return uniprot_entry
            elif action == "start" and tag == self.tag_accession:
                prev_text = elem.text
                if elem.text:
                    uniprot_entry["accession"].append(elem.text)
            elif action == "end" and tag == self.tag_accession:
                if prev_text is None and elem.text:
                    uniprot_entry["accession"].append(elem.text)
            elif action == "start" and tag == self.tag_name:
                uniprot_entry["name"].append(elem.text)
            elif action == "end" and tag == self.tag_name:
                pass
            else:
                print("no capture", action, tag, elem.attrib, elem.text)
