#!/usr/bin/env bash

for line in $(ls -1 gatk4_json/picard_vcf_*.json); do
    python3 parse_gatk_json.py --json ${line} --json_type picard_vcf --xml_out /Users/letaw/PycharmProjects/gatk4_xml_prepare/output
done
