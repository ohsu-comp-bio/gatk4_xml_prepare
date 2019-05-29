#!/usr/bin/env python

from lxml import etree

handle_out = open('blahblah', 'w')

name='gatk4_vcftointervallist'

tool = etree.Element("tool", id=name, name="GATK4 VcfToIntervalList (Picard)", version="@WRAPPER_VERSION@0", profile="17.09")
tool.append( etree.Element("child1") )
handle_out.write(etree.tostring(tool, pretty_print=True, encoding = "unicode"))

handle_out.close()