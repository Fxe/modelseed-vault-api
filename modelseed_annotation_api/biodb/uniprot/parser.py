import logging
import lxml.etree as ET

logging.getLogger(__name__)


class UnirefParser:
    def __init__(self, parse_limit=None):
        """

        :param parse_limit: max number of entries
        """
        self.max_entries = parse_limit

        self.tag_entry = "entry"
        self.tag_property = "property"
        self.tag_sequence = "sequence"
        self.tag_db_reference = "dbReference"
        self.tag_representative_member = "representativeMember"
        self.tag_member = "member"
        self.tag_name = "name"

    def set_xml_ns(self, ns):
        self.tag_entry = ns + self.tag_entry
        self.tag_property = ns + self.tag_property
        self.tag_sequence = ns + self.tag_sequence
        self.tag_db_reference = ns + self.tag_db_reference
        self.tag_member = ns + self.tag_member
        self.tag_representative_member = ns + self.tag_representative_member
        self.tag_name = ns + self.tag_name

    def parse(self, fh):
        yielded = 0
        parser = ET.iterparse(fh, events=("end", "start"))
        for action, elem in parser:
            tag = elem.tag
            # print(action, tag)
            if action == "start" and tag == self.tag_entry:
                yield self.parse_entry(elem, parser)
                yielded += 1
            else:
                print("parse", action, elem)
            if self.max_entries and yielded >= self.max_entries:
                break

    def parse_sequence(self, elem, parser):
        sequence = dict(elem.attrib)
        sequence["value"] = elem.text
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.tag_sequence:
                return sequence
            elif action == "start" and tag == self.tag_property:
                k = elem.attrib["type"]
                value = elem.attrib["value"]
                if k not in sequence:
                    sequence[k] = []
                sequence[k].append(value)
            else:
                print("parse_sequence", action, elem)
        print("ERROR!")

    def parse_db_reference(self, elem, parser):
        db_reference = dict(elem.attrib)
        for action, elem in parser:
            tag = elem.tag
            if action == "end" and tag == self.tag_db_reference:
                return db_reference
            elif action == "start" and tag == self.tag_property:
                k = elem.attrib["type"]
                value = elem.attrib["value"]
                if k not in db_reference:
                    db_reference[k] = []
                db_reference[k].append(value)
            elif action == "end" and tag == self.tag_property:
                pass
            else:
                print("parse_db_reference", action, elem)
        print("ERROR!", elem.sourceline)

    def parse_member(self, elem, parser, end_tag):
        uniref_entry_member = {"db_reference": []}
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and tag == self.tag_sequence:
                sequence = self.parse_sequence(elem, parser)
                if "sequence" not in uniref_entry_member:
                    uniref_entry_member["sequence"] = sequence
                else:
                    print("double sequence", elem.sourceline)
            elif action == "start" and tag == self.tag_db_reference:
                db_reference = self.parse_db_reference(elem, parser)
                uniref_entry_member["db_reference"].append(db_reference)
            elif action == "start" and tag == self.tag_property:
                k = elem.attrib["type"]
                value = elem.attrib["value"]
                if k not in uniref_entry_member:
                    uniref_entry_member[k] = []
                uniref_entry_member[k].append(value)
            elif action == "end" and tag == self.tag_property:
                pass
            elif action == "end" and tag == end_tag:
                return uniref_entry_member
            else:
                print("parse_representative_member", action, elem)
        print("ERROR!", elem.sourceline)

    def parse_entry(self, elem, parser):
        uniref_entry = dict(elem.attrib)
        uniref_entry["representative_members"] = []
        uniref_entry["members"] = []
        for action, elem in parser:
            tag = elem.tag
            if action == "start" and tag == self.tag_name:
                uniref_entry["name"] = elem.text
            elif action == "end" and tag == self.tag_name:
                pass
            elif action == "start" and tag == self.tag_property:
                k = elem.attrib["type"]
                value = elem.attrib["value"]
                if k not in uniref_entry:
                    uniref_entry[k] = []
                uniref_entry[k].append(value)
            elif action == "end" and tag == self.tag_property:
                pass
            elif action == "start" and tag == self.tag_representative_member:
                uniref_entry_representative_member = self.parse_member(
                    elem, parser, self.tag_representative_member
                )
                uniref_entry["representative_members"].append(
                    uniref_entry_representative_member
                )
            elif action == "start" and tag == self.tag_member:
                uniref_entry_member = self.parse_member(elem, parser, self.tag_member)
                uniref_entry["members"].append(uniref_entry_member)
            elif action == "end" and tag == self.tag_entry:
                return uniref_entry
            else:
                print("parse_entry", action, elem)
