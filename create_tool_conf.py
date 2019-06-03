# Create a tool_conf file based on all jsons in a particular directory.

from lxml import etree
from string import Template
import argparse
import json
import os

VERSION = '0.0.1'

def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--json_dir', help='Input JSON Dir')
    parser.add_argument('--tool_conf', help='Output tool_conf.xml')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

class XmlFile(object):

    def __init__(self, filename, prefix='gatk4_auto_dev'):
        self.filename = self._filename_create(filename)
        self.prefix = prefix
        self.xml_tmpl = Template('${prefix}/${prefix}_${fileid}.xml')
        self.full_path = self.xml_tmpl.substitute(prefix=self.prefix, fileid=self.filename)

    def _filename_create(self, filename):
        """
        Return the name, but reformatted.
        :param prefix:
        :return:
        """
        return filename.lower().split(' ')[0]


class ToolConfWrite():
    """
    Write, using etrees.
    """
    def __init__(self, filenames, output, sid='gatk4_auto_dev', stitle='GATK4 Auto Generated Tools'):
        self.stitle = stitle
        self.sid = sid
        self.output = output
        self.filenames = filenames
        self.to_write = etree.tostring(self._assemble_etree(), pretty_print=True, encoding="unicode")

    def _assemble_etree(self):
        """
        <section id="fgbio" name="fgbio">
        <tool file="gatk4_auto_dev/gatk4_auto_dev_annotatepairorientation.xml"/>
        :return:
        """
        section = etree.Element('section', id=self.sid, name=self.stitle)
        for filename in self.filenames:
            etree.SubElement(section, 'tool', file=filename)
        return section

    def write_me(self):
        """
        Write to file.
        :return:
        """
        handle_out = open(self.output, 'w')
        handle_out.write(self.to_write)
        handle_out.close()


def main():
    args = supply_args()
    all_files = []

    print(os.listdir(args.json_dir))
    for entry in os.scandir(args.json_dir):
        with open(entry, 'rU') as myjson:
            json_file = json.load(myjson)
            toolname = json_file['name']
            all_files.append(XmlFile(toolname).full_path)

    ToolConfWrite(all_files, args.tool_conf).write_me()


if __name__ == "__main__":
    main()
