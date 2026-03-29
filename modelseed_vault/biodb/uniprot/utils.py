import logging

logger = logging.getLogger(__name__)

HEADER = """<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/docs/uniprot.xsd">
"""
FOOTER = "</uniprot>"


def split_into_blocks(
    fh_input, output_folder, l_block_size=2500, l_max=None, h=HEADER, f=FOOTER, decode=None
):
    i = 0
    line = fh_input.readline()
    if decode:
        line = line.decode(decode)
    xml_record = None
    block_index = 0
    written = 0
    fh_write = open(f"{output_folder}/block_{block_index}_{l_block_size}.xml", "w")
    if h:
        fh_write.write(h)
    while line:
        line = fh_input.readline()
        if decode:
            line = line.decode(decode)

        i += 1

        if "<entry dataset" in line and xml_record is None:
            xml_record = ""
            xml_record += line
        elif "</entry" in line:
            xml_record += line
            fh_write.write(xml_record)
            xml_record = None
            written += 1
            if written >= l_block_size:
                if f:
                    fh_write.write(f)
                fh_write.close()
                block_index += 1
                written = 0
                fh_write = open(
                    f"{output_folder}/block_{block_index}_{l_block_size}.xml", "w"
                )
                if h:
                    fh_write.write(h)
        else:
            if xml_record:
                xml_record += line

        if l_max and i > l_max:
            logger.warning("max lines exceeded term")

    if f:
        fh_write.write(f)
    fh_write.close()
    return block_index
