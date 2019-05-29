#!/usr/bin/env python

# Potentially make List[String] types a repeat xml element.  Otherwise, it appears these are separated by commas on the command line.
# Annotation options need macros.
# Macros for read filter options?  These have their own json files it seems...
# gVCF not currently set to be different than VCF.

from lxml import etree
from string import Template
from xml.sax.saxutils import escape
import argparse
import json
import pypandoc

VERSION="0.3.0"

def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--json', help='Input JSON')
    parser.add_argument('--xml_out', help='Output Directory')
    # parser.add_argument('--old_galaxy', action="store_true", help="Produce XML for Galaxy versions that don't support section tag")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class Mappings(object):
    """
    Hold all of the mappings we need to provide for parameters.
    """
    def __init__(self):
        # These are parameters that are already covered in the macros.xml file.  They also may be options that we
        # don't want included in any GATK wrapper.
        # REFERENCE_SEQUENCE: If reference is denoted optional then it should go in the optional section.  Build function to place it.
        # Should be able to utilize self.blob['required'] to determine section (either optional or main), then look at self.macro_to_chth to determine
        # the macro that REFERENCE_SEQUENCE is associated to.
        self.xml_json_type_map = {
                             'Boolean': 'boolean',
                             'boolean': 'boolean',
                             'Long': 'integer',
                             'long': 'integer',
                             'Integer': 'integer',
                             'int': 'integer',
                             'byte': 'integer',
                             'Double': 'float',
                             'double': 'float',
                             'Float': 'float',
                             'float': 'float',
                             'Format': 'select',
                             'String': 'text',
                             'Set[String]': 'text',
                             'SortOrder': 'select',
                             'File': 'data',
                             'Set[File]': 'data',
                             'Set[MetricAccumulationLevel]': 'text',
                             'List[File]': 'data',
                             'List[Double]': 'text',
                             'List[Integer]': 'integer',
                             'List[String]': 'text',
                             'List[Type]': 'text',
                             'List[FeatureInput[VariantContext]]': 'data',
                             'ArrayList[String]': 'text',
                             'FeatureInput[VariantContext]': 'data',
                             'LogLevel': 'select',
                             'Implementation': 'select',
                             'IntervalMergingRule': 'select',
                             'IntervalSetRule': 'select',
                             'ReferenceConfidenceMode': 'select',
                             'ValidationStringency': 'select',
                             'GenotypingOutputMode': 'select',
                             'OutputMode': 'select',
                             'PCRErrorModel': 'select',
                             'WriterType': 'select',
                             'NumberAlleleRestriction': 'select'
                             }
        self.xml_json_num_map = {'Infinity': '', '-Infinity': '', 'NA': ''}
        self.xml_json_req_map = {'no': 'true', 'yes': 'false'}
        # If the parameter is in here, it won't be included in either the cheetah section or the param section.
        # self.common_args = ("version", "showHidden", "help", "arguments_file", "VERBOSITY",
        #                     "VALIDATION_STRINGENCY", "USE_JDK_INFLATER", "USE_JDK_DEFLATER", "TMP_DIR", "QUIET",
        #                     "MAX_RECORDS_IN_RAM", "GA4GH_CLIENT_SECRETS", "CREATE_MD5_FILE", "CREATE_INDEX",
        #                     "COMPRESSION_LEVEL", "REFERENCE_SEQUENCE", "OUTPUT", "SEQUENCE_DICTIONARY", "INPUT",
        #                     "input", "reference", "output", "annotation", "annotation_group", "annotations_to_exclude",
        #                     "intervals", "exclude_intervals", "read_index", "interval_padding",
        #                     "interval_exclusion_padding", "output_prefix", "sequence_dictionary", "variant")
        self.read_filter_args = ('')
        # There seems to be a significant amount of information that can't be retrieved from the provided json.
        # Make this in to a config file.
        # Might be able to classify a bunch of this based on input output files for a particular tool.  For instance, vcf_tabix
        # is applicable if the input is a VCF.
        # Optional status of a partiular parameters could be connected to macros as well.  For instance, ref_sel is the macro for the reference genome.
        # So, the 'ref_sel': 'reference' relation could be used along with 'required' status to place these.
        # One challenge with this is the macros would have to either be prebuilt with the necessary section variable, or
        # would have to be rewritten on the fly.

        # self.tool_data = {'SortVcf':
        #                       {'output_fmt': {},
        #                        'pre_tmpls': ['vcf_tabix_multi'],
        #                        'post_tmpls': ['picard_opts', 'picard_vcf_output_opts', 'vcf_input_multi', 'picard_ref_opts_opt', 'picard_seqdict_opts'],
        #                        'pre_params': ['gzip_vcf_params', 'vcf_input_params_multi'],
        #                        'opt_params': ['ref_sel', 'seq_dict_sel'],
        #                        'post_params': ['picard_params'],
        #                        'output_params': ['gzip_vcf_output_params', 'picard_output_params']},
        #                   'GenotypeGVCFs':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': ['gatk_tabix', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
        #                        'post_tmpls': ['ref_opts', 'vcf_output_opts', 'gatk_input', 'gatk_ints_chth', 'gatk_excl_ints_chth'],
        #                        'pre_params': ['ref_sel', 'gzip_vcf_params', 'gatk_ints', 'gatk_vcf_input_params',
        #                                       'gatk_excl_ints'],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['gzip_vcf_output_params']},
        #                   'CombineGVCFs':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': ['gatk_tabix_multi', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
        #                        'post_tmpls': ['ref_opts', 'vcf_output_opts', 'gatk_input_multi', 'gatk_ints_chth', 'gatk_excl_ints_chth'],
        #                        'pre_params': ['ref_sel', 'gzip_vcf_params', 'gatk_vcf_input_params_multi', 'gatk_ints', 'gatk_excl_ints'],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['gzip_vcf_output_params']},
        #                   'HaplotypeCaller':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': ['bam_index', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
        #                        'post_tmpls': ['ref_opts', 'vcf_output_opts', 'gatk_bam_input', 'gatk_ints_chth', 'gatk_excl_ints_chth'],
        #                        'pre_params': ['ref_sel', 'gatk_req_params', 'gzip_vcf_params', 'gatk_ints', 'gatk_excl_ints'],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['gzip_vcf_output_params']},
        #                   'SelectVariants':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': ['gatk_tabix', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
        #                        'post_tmpls': ['ref_opts', 'vcf_output_opts', 'gatk_input', 'gatk_ints_chth', 'gatk_excl_ints_chth'],
        #                        'pre_params': ['ref_sel', 'gzip_vcf_params', 'gatk_vcf_input_params', 'gatk_ints', 'gatk_excl_ints'],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['gzip_vcf_output_params']},
        #                   'CollectReadCounts':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': ['bam_index', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
        #                        'post_tmpls': ['ref_opts_opt', 'gatk_bam_input', 'gatk_ints_chth', 'gatk_excl_ints_chth', 'hdf5_output_chth'],
        #                        'pre_params': ['gatk_req_params', 'gatk_ints', 'gatk_excl_ints'],
        #                        'opt_params': ['ref_sel'],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['hdf5_output']},
        #                   'PlotDenoisedCopyRatios':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': ['bam_index', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
        #                        'post_tmpls': ['ref_opts_opt', 'vcf_output_opts', 'gatk_bam_input', 'gatk_ints_chth', 'gatk_excl_ints_chth'],
        #                        'pre_params': ['ref_sel', 'gatk_req_params', 'gzip_vcf_params', 'gatk_ints',
        #                                       'gatk_excl_ints'],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['gzip_vcf_output_params']},
        #                   'ModelSegments':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': [],
        #                        'post_tmpls': ['modelsegments_chth'],
        #                        'pre_params': [],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['modelsegments_output']},
        #                   'CallCopyRatioSegments':
        #                       {'output_fmt': {'output': 'tabular'},
        #                        'input_fmt': {'input': 'tabular'},
        #                        'pre_tmpls': [],
        #                        'post_tmpls': [],
        #                        'pre_params': [],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': []},
        #                   'PlotModeledSegments':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': [],
        #                        'post_tmpls': ['plotmodeledsegments_chth', 'gatk_seqdict'],
        #                        'pre_params': ['seq_dict_sel'],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['plotmodeledsegments_output']},
        #                   'CollectAllelicCounts':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': ['bam_index', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
        #                        'post_tmpls': ['ref_opts_opt', 'vcf_output_opts', 'gatk_bam_input', 'gatk_ints_chth', 'gatk_excl_ints_chth'],
        #                        'pre_params': ['ref_sel', 'gatk_req_params', 'gzip_vcf_params', 'gatk_ints', 'gatk_excl_ints'],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': ['gzip_vcf_output_params']},
        #                   'DenoiseReadCounts':
        #                       {'output_fmt': {'denoised_copy_ratios': 'tabular',
        #                                       'standardized_copy_ratios': 'tabular'},
        #                        'input_fmt': {},
        #                        'pre_tmpls': [],
        #                        'post_tmpls': ['hdf5_input_chth'],
        #                        'pre_params': ['hdf5_input'],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': []},
        #                   'CollectHsMetrics':
        #                       {'output_fmt': {'OUTPUT': 'tabular'},
        #                        'input_fmt': {'INPUT': 'sam,bam'},
        #                        'pre_tmpls': ['picard_bam_index'],
        #                        'post_tmpls': ['picard_opts'],
        #                        'pre_params': [],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': ['picard_params'],
        #                        'output_params': []},
        #                   'Mutect2':
        #                       {'output_fmt': {},
        #                        'input_fmt': {},
        #                        'pre_tmpls': [],
        #                        'post_tmpls': [],
        #                        'pre_params': [],
        #                        'opt_params': [],
        #                        'adv_params': [],
        #                        'post_params': [],
        #                        'output_params': []}
        #                   }


        # For a parameter to be included here, it must have a consistent file format throughout all GATK4 JSONs.
        # Otherwise, this will need to be specified as a macro, and macro will need to be connected to tool.
        self.gen_out_fmt = {'activity_profile_out': 'tabular',
                            'assembly_region_out': 'tabular',
                            'graph_output': 'txt',
                            'bam_output': 'bam',
                            }

        # This is meant to be for those parameters that are booleans, but instruct for an output file to be created.
        # This is a work in progress, will probably ignore these parameters for the time being.
        self.out_create_params = {}

        self.gen_in_fmt = {'alleles': 'vcf',
                           'pedigree': 'tabular',
                           'population_callset': 'vcf',
                           'panel_of_normals': 'vcf',
                           'contamination_fraction_per_sample_file': 'tabular',
                           'concordance': 'vcf',
                           'comp': 'vcf',
                           'discordance': 'vcf',
                           'gatk_config_file': 'txt',
                           'germline_resource': 'vcf',
                           'annotated_intervals': 'gatk_interval',
                           'count_panel_of_normals': 'h5',
                           'allelic_counts': 'tabular',
                           'normal_allelic_counts': 'tabular',
                           'denoised_copy_ratios': 'tabular',
                           'dbsnp': 'vcf',
                           'BAIT_INTERVALS': 'picard_interval_list',
                           'TARGET_INTERVALS': 'picard_interval_list',
                           'segments': 'tabular',
                           'intervals': 'gatk_interval',
                           'read_index': 'idx'
                           }

        # This is for the parameters in the XML
        self.param_to_macro_xml = {'exclude_intervals': ['gatk_excl_ints'],
                                   'intervals': ['gatk_ints'],
                                   'REFERENCE_SEQUENCE': ['ref_sel'],
                                   'REFERENCE': ['ref_sel'],
                                   'reference': ['ref_sel'],
                                   'SEQUENCE_DICTIONARY': ['seq_dict_sel'],
                                   'sequence_dictionary': ['seq_dict_sel']
                                   }

        # Template {<pname>: ((<main_chth>), (<pre_chth>))}
        self.param_to_macro_tmpl = {'exclude_intervals': ['gatk_excl_ints_chth'],
                                    'intervals': ['gatk_ints_chth'],
                                    'REFERENCE_SEQUENCE': ['picard_ref_opts'],
                                    'REFERENCE': ['picard_ref_opts_plain'],
                                    'reference': ['ref_opts'],
                                    'SEQUENCE_DICTIONARY': ['picard_seqdict_opts'],
                                    'sequence_dictionary': ['gatk_seqdict']
                                    }

        self.param_to_macro_pretmpl = {'exclude_intervals': ['pre_gatk_excl_ints_chth'],
                                       'intervals': ['pre_gatk_ints_chth']
                                       }

        # self.param_to_macro_tmpl = {'REFERENCE_SEQUENCE': {'optional': ['picard_ref_opts_opt'], 'required': ['picard_ref_opts']},
        #                             'REFERENCE': {'required': ['picard_ref_opts_plain']},
        #                             'reference': {'optional': ['ref_opts_opt'], 'required': ['ref_opts']}
        #                             }

        self.param_tmpls = {'integer': ['name', 'argument', 'type', 'optional', 'value', 'min', 'max', 'label', 'help'],
                            'float': ['name', 'argument', 'type', 'optional', 'value', 'min', 'max', 'label', 'help'],
                            'text': ['name', 'argument', 'type', 'optional', 'value', 'label', 'help'],
                            'data': ['name', 'argument', 'type', 'optional', 'format', 'label', 'help'],
                            'select': ['name', 'argument', 'type', 'optional', 'label', 'help'],
                            'boolean': ['name', 'argument', 'type', 'truevalue', 'falsevalue', 'optional', 'checked', 'label', 'help'],
                            'output': ['format', 'name', 'label', 'help']}

        # Many parameters will take multiple Galaxy datatypes.  Providing a mapping will give flexibility to add or
        # subtract types when wanted.
        self.file_type_map = {'vcf': 'vcf,vcf_bgzip',
                              'idx': 'tabix,bai',
                              'gatk_interval': 'gatk_interval,bed,vcf',
                              'sam': 'sam,bam',
                              'tabular': 'tabular'}

        # Relate tool name to a common argument name, such as input (which can be connected to different file types),
        # and a file type.
        self.tool_file_type = {'Mutect2': {'input': 'sam'},
                               'SortSam': {'INPUT': 'sam'}}
        self.tool_output_file_type = {'Mutect2': {'output': 'vcf'},
                                      'SortSam': {'OUTPUT': 'sam'}}

        # Relate the name of the argument to the selection of macros that goes with it.
        self.macro_to_param = {'input': {'sam': {'pre_chth': ['bam_index_pre_chth'],
                                                 'main_chth': ['gatk_bam_req_chth'],
                                                 'main_xml': ['gatk_bam_req_params']}},
                               'INPUT': {'sam': {'pre_chth': ['picard_bam_index'],
                                                 'main_chth': ['picard_bam_input'],
                                                 'main_xml': ['gatk_bam_req_params']}},
                               'output': {'vcf': {'main_chth': ['vcf_output_opts'],
                                                  'main_xml': ['gzip_vcf_params'],
                                                  'out_xml': ['gzip_vcf_output_params']}},
                               'OUTPUT': {'sam': {'main_chth': [],
                                                  'main_xml': [],
                                                  'out_xml': []}}
                               }

        # For parameters we universally do not want to see in the UI.
        # interval_padding is taken care of any time we see an intervals parameter.
        self.ignore_params = ('help', 'version', 'showHidden', 'interval_padding', 'interval_exclusion_padding',
                              'TMP_DIR',
                              'create_output_bam_index', 'create_output_bam_md5', 'read_index',
                              'create_output_variant_index', 'create_output_variant_md5')

    def map_lookup(self, val, dict):
        """

        :param self:
        :return:
        """
        try:
            return dict[val]
        except:
            return val


class XmlTemplates(object):
    def __init__(self):
        # Cheetah section
        self.chth_tmpl = PercentTemplate('#include source=$%macro#')
        # Required argument.
        self.req_out_chth = PercentTemplate('%argument $%name')
        # Required boolean argument.
        self.req_out_chth_bool = PercentTemplate('$%name')
        self.vcf_choose = PercentTemplate('#if $%section.%name'
                                          '\n#if $%section.%name.is_of_type("vcf_bgzip")'
                                          '\n%argument %name.vcf.gz'
                                          '\n#else'
                                          '\n%argument %name.vcf'
                                          '\n#end if'
                                          '\n#end if')
        self.vcf_tabix = PercentTemplate('#if $%section.%name'
                                         '\n#set datatype = $%section.%name.datatype'
                                         '\n#if $%section.%name.is_of_type("vcf_bgzip")'
                                         '\nln -s $%section.%name %name.vcf.gz &&'
                                         '\ntabix %name.vcf.gz &&'
                                         '\n#else'
                                         '\nln -s $%section.%name %name.vcf &&'
                                         '\n#end if'
                                         '\n#end if')
        self.file_chth = PercentTemplate('#if $%section.%out_sel_name\n%argument $%name\n#end if')
#        self.file_chth_old_gal = PercentTemplate('#if str($output_opt.output_opt_sel) == "yes":\n#if $output_opt.%out_sel_name:\n%argument $%name\n#end if\n#end if')
        self.ext_arg = PercentTemplate('#if $%section.%name\n  %argument $%section.%name\n#end if\n')
        self.ext_arg_bool = PercentTemplate('#if $%section.%name\n  $%section.%name\n#end if\n')
#        self.ext_arg_old_gal = '#if str($%section.%{section}_sel) == "yes":\n#if $%section.%name:\n%argument $%section.%name\n#end if\n#end if'
        self.reg_arg = '#if $%name\n  %argument $%name\n#end if\n'

        # XML section. Most of this is handled via etrees, but some cases are more easily handled here.
        self.xml_tmpl = Template('<expand macro="$macro_name" />')
        self.sel_out = Template('<option value="$value" selected="$selected">$value</option>')
        self.out_label_tmpl = Template('$${tool.name} on $${on_string}: $name $format')
        self.out_tmpl = Template('<data format="$format" name="$name" label="$${tool.name} on $${on_string}: $format" />')


class CheetahPrep(XmlTemplates):
    """
    Take arguments, then prep the Cheetah string.
    Rules:
    is_req (req_out_chth) (no section)
    is_req && has_macro (chth_tmpl) (no section)
    is_req && has_pre_macro (chth_tmpl) (no section)

    all non-req args go in sections

    not_req (ext_arg)
    not_req && has_macro (chth_tmpl) (refer to section macro)
    not_req && has_pre_macro (chth_tmpl) (refer to section macro)

    Need pname, macro name

    """
    def __init__(self, pname, argname, section, is_req=False, mname=None, pre_mname=None, is_bool=False, out_sel_name=None):
        XmlTemplates.__init__(self)
        # Output is required, and there is a macro in my_xml.macros_tmpl.

        self.pname = pname
        self.argname = argname
        self.section = section
        self.mname = mname
        self.pre_mname = pre_mname
        self.is_req = is_req
        self.is_bool = is_bool
        self.out_sel_name = out_sel_name

        self.chth = self._chth_create()

    def _create_section_macro(self, macro):
        """
        Create a macro name appropriate for macros listed in sections.
        :param macro:
        :return:
        """
        return '_'.join([macro, self.section])

    def _chth_create(self):
        """
        Logic to decide which Cheetah template to return.
        :return:
        """
        if self.is_req and self.mname and not self.pre_mname:
            return self.chth_tmpl.substitute(macro=self.mname)

        if self.is_req and not self.mname and not self.pre_mname and self.is_bool:
            return self.req_out_chth_bool.substitute(name=self.pname)

        if self.is_req and not self.mname and not self.pre_mname and not self.is_bool:
            return self.req_out_chth.substitute(argument=self.argname, name=self.pname)

        if self.is_req and not self.mname and self.pre_mname:
            return self.chth_tmpl.substitute(macro=self.pre_mname)

        if not self.is_req and self.out_sel_name:
            return self.file_chth.substitute(out_sel_name=self.out_sel_name, section=self.section,
                                             name=self.pname, argument=self.argname)

        if not self.is_req and not self.mname and not self.pre_mname and self.is_bool:
            return self.ext_arg_bool.substitute(section=self.section, name=self.pname)

        if not self.is_req and not self.mname and not self.pre_mname and not self.is_bool:
            return self.ext_arg.substitute(section=self.section, argument=self.argname, name=self.pname)

        if not self.is_req and self.mname and not self.pre_mname:
            # return self.chth_tmpl.substitute(macro=self._create_section_macro(self.mname))
            return self.chth_tmpl.substitute(macro=self.mname)

        if not self.is_req and not self.mname and self.pre_mname:
            return self.chth_tmpl.substitute(macro=self.pre_mname)

        if self.mname and self.pre_mname:
            raise Exception("You can't pass both mname and pre_mname to the same CheetahPrep object.")


class JsonXml(Mappings, XmlTemplates):
    """
      This class will take one json-formatted entry from the "arguments" section of a GATK4 json file, and will provide
      the variables we will need to do the rest of the stuff we have to do when preparing an XML file.

      Under arguments, dict's look like this:

      "summary": "read one or more arguments files and add them to the command line",
      "name": "--arguments_file",
      "synonyms": "NA",
      "type": "List[File]",
      "required": "no",
      "fulltext": "",
      "defaultValue": "[]",
      "minValue": "NA",
      "maxValue": "NA",
      "minRecValue": "NA",
      "maxRecValue": "NA",
      "kind": "optional",
      "options": []

      if required then this param should go in the main portion of the xml, not in a section.
      macros.xml entries based off of (name, type, kind)?

      if REFERENCE_SEQUENCE or REFERENCE or reference, attach appropariate macro based on section designation
      required will not go in a section
      make sure data type is "data"
      also, "String" types will need the same macros for reference sequence, only see String and File
      can be common, required, optional, so we need macros for each of these, unless there is a better way to substitute section
      on the fly
      could also force both optional and common in to the same optional section, might make it easier for the user

    """
    def __init__(self, blob, tool_name, args):
        """
        :param blob:
        blob == {'summary': 'display the version number for this tool', 'name': '--version', 'synonyms': 'NA', 'type': 'boolean', 'required': 'no', 'fulltext': '', 'defaultValue': 'false', 'minValue': 'NA', 'maxValue': 'NA', 'minRecValue': 'NA', 'maxRecValue': 'NA', 'kind': 'optional', 'options': []}
        """
        # Access mappings and templates for each parameter.
        Mappings.__init__(self)
        XmlTemplates.__init__(self)
        # This comes from the json file provided from GATK4.
        self.args = args
        self.blob = blob
        # Tool name, currently passed as input.
        self.tool_name = tool_name
        # Type comes from mapping of json 'types' to types recognized in a Galaxy xml.
        try:
            self.type = self.xml_json_type_map[blob['type']]
        except:
            raise Exception('Argument type %s not recognized.' % self.blob['type'])
        # Make parameter name available more readily.
        self.pname = self.blob['name'].lstrip('-').replace('-', '_')
        # The 'kind' category determines section, and required status of parameter.
        # Potential values are ('advanced', 'common', 'deprecated', 'optional', 'positional', 'required')
        self.section = blob['kind']
        # If argument kind value is defined as common, set self.common as True.
        self.common = (self.section == 'common')
        # Since output status is not listed in the json blob, we provide it as a mapping.
        # TODO: Potentially search for phrase "Output" or "output" from description and set status accordingly.
        self.is_output = ((self.pname in self.gen_out_fmt) or
                          (self.pname in self.tool_output_file_type[self.tool_name]) or
                          (self.pname in self.out_create_params))
        # Define this parameter as required output.
        if self.is_output and self.section == 'required':
            self.is_req_output = True
        else:
            self.is_req_output = False
        # If this is not required output, define additional variables to allow for selecting whether or not we
        # will see a particular output.
        if self.is_output and not self.is_req_output and self.pname in self.gen_out_fmt:
            self.section = 'output_opt'
            self.out_sel_name = self.pname + '_sel'
            self.out_sel_arg = self.blob['name'] + '_sel'
        elif self.is_output and not self.is_req_output and self.pname in self.out_create_params:
            self.section = 'output_opt'
            self.out_sel_name = self.pname
            self.out_sel_arg = self.blob['name']
        else:
            self.out_sel_name = None
            self.out_sel_arg = None

        # If this parameter is defined as 'data' from our mapping, and was not defined as output, it is an input file.
        self.is_input = (not self.is_output) and self.type == 'data' \
                        or self.pname in self.tool_file_type[self.tool_name] \
                        or self.pname in self.gen_in_fmt
        if self.is_input and self.section == 'required':
            self.is_req_input = True
        else:
            self.is_req_input = False

        if self.section == 'required':
            self.is_req = True
        else:
            self.is_req = False


        # This will hold a new structure, corresponding to xml key/value pairs.
        self.xml_out = self.reblob()

        if self.xml_out['type'] == 'boolean':
            self.is_bool = True
        else:
            self.is_bool = False

        # if self.pname in self.gen_in_fmt:
        #     self.is_input_vcf = self.is_input and self.gen_in_fmt[self.pname] == 'vcf,vcf_bgzip'
        # else:
        #     self.is_input_vcf = None

        # Set to common if this argument is seen inside the known common arguments.
        # self.common = self.pname not in self.tool_data[self.tool_name]['output_fmt'] and \
        #               self.pname not in self.tool_data[self.tool_name]['input_fmt']
        #self.has_mcro_xml = self.xml_out['name'] in self.macro_xml
        #self.has_mcro_tmpl = self.xml_out['name'] in self.macro_tmpl
        # Check to see if the type, as defined in the json file, is recognized.
        # if self.blob['type'] not in self.xml_json_type_map:
        #     print('Argument type %s not recognized, skipping.' % self.blob['type'])

        # If there are options listed for an argument, set this as a select argument, which will need additional handling
        # when writing to XML.
        if self.blob['options']:
            self.is_select = True
            self.sel_blob = self.sel_prep()
        else:
            self.is_select = False
            self.sel_blob = None

        self.is_input_vcf = None

        self.macros = {'main_chth': [],
                       'pre_chth': [],
                       'main_xml': [],
                       'out_xml': []}

        # self.chth_pre = self.macros['pre_chth']

        # Set input or output types, and provide dictionary of macro types for this parameter/tool.
        if self.tool_name in self.tool_file_type:
            if self.pname in self.tool_file_type[self.tool_name]:
                self.input_type = self.tool_file_type[self.tool_name][self.pname]
                self._macro_update(self.macro_to_param[self.pname][self.input_type])
            elif self.pname in self.tool_output_file_type[self.tool_name]:
                self.output_type = self.tool_output_file_type[self.tool_name][self.pname]
                self._macro_update(self.macro_to_param[self.pname][self.output_type])

        if self.pname in self.param_to_macro_xml:
            macros_xml = self.param_to_macro_xml[self.pname]
        elif self.pname in self.macro_to_param and self.is_input:
            self.raw_format = self.tool_file_type[self.tool_name][self.pname]
            macros_xml = self.macro_to_param[self.pname][self.raw_format]["main_xml"]
        elif self.pname in self.macro_to_param and self.is_output:
            self.raw_format = self.tool_output_file_type[self.tool_name][self.pname]
            macros_xml = self.macro_to_param[self.pname][self.raw_format]["main_xml"]
        else:
            macros_xml = None
        if macros_xml:
            self.macros['main_xml'].extend(macros_xml)

        if self.pname in self.param_to_macro_tmpl:
            self._simple_macro_update(self.param_to_macro_tmpl[self.pname], 'main_chth')
            # self.macros['main_chth'].extend(self.param_to_macro_tmpl[self.pname][self.section])

        if self.pname in self.param_to_macro_pretmpl:
            self._simple_macro_update(self.param_to_macro_pretmpl[self.pname], 'pre_chth')

        # if self.pname in self.macro_to_param:
        #     self.main_chth = None

        # Create Cheetah strings.
        self.chth = []
        if self.macros['main_chth']:
            for macro in self.macros['main_chth']:
                self.chth.append(CheetahPrep(self.pname, self.xml_out['argument'], self.section, self.is_req, mname=macro, is_bool=self.is_bool, out_sel_name=self.out_sel_name).chth)
        else:
            self.chth.append(CheetahPrep(self.pname, self.xml_out['argument'], self.section, self.is_req, mname=None,
                                         is_bool=self.is_bool, out_sel_name=self.out_sel_name).chth)

        self.chth_pre = []
        if self.macros['pre_chth']:
            for macro in self.macros['pre_chth']:
                self.chth_pre.append(CheetahPrep(self.pname, self.xml_out['argument'], self.section, self.is_req, pre_mname=macro, is_bool=self.is_bool).chth)

        # Fill in XML tags to dict.
        if self.is_output:
            self.xml_param = {label:self.xml_out[label] for label in self.param_tmpls['output']}
        else:
            if not self.macros['main_xml']:
                self.xml_param = {label:self.xml_out[label] for label in self.param_tmpls[self.xml_out['type']]}
            else:
                self.xml_param = None

        # if self.pname in self.param_to_macro_tmpl:
        #     self.macros_tmpl = self.param_to_macro_tmpl[self.pname][self.section]
        #     self.chth = self.chth_tmpl.substitute(macro=self.macros_tmpl)
        # else:
        #     self.macros_tmpl = None

        # if self.pname in self.param_to_macro_pretmpl:
        #     self.macros_pretmpl = self.param_to_macro_pretmpl[self.pname][self.section]
        # else:
        #     self.macros_pretmpl = None

        # if self.pname == 'create_output_bam_index':
        # self._params_stdout()
        # print('\n')

    def _params_stdout(self):
        """
        Print variable contents to stdout.
        :return:
        """
        print("JSON BLOB: {0}".format(self.blob))
        print("Tool Name: {0}".format(self.tool_name))
        print("Type: {0}".format(self.type))
        print("Parameter Name: {0}".format(self.pname))
        print("XML Section: {0}".format(self.section))
        print("XML Out: {0}".format(self.xml_out))

        print("Output Select Name: {0}".format(self.out_sel_name))
        print("Output Select Arg: {0}".format(self.out_sel_arg))
        print("Select Blob: {0}".format(self.sel_blob))

        print("Cheetah String: {0}".format(self.chth))
        print("Pre Cheetah String: {0}".format(self.chth_pre))
        # print("Cheetah Macros: {0}".format(self.macros_tmpl))
        # print("Cheetah Pre Macros: {0}".format(self.macros_pretmpl))

        print("XML Params: {0}".format(self.xml_param))
        print("Macros: {0}".format(self.macros))

        print("Is Common: {0}".format(self.common))
        print("Is Input: {0}".format(self.is_input))
        print("Is Req Input: {0}".format(self.is_req_input))
        print("Is Output: {0}".format(self.is_output))
        print("Is Req Output: {0}".format(self.is_req_output))
        print("Is Select: {0}".format(self.is_select))

    def _macro_update(self, new_macros):
        """
        Update the macros structure.  Structure should look like:
        self.macros = {'main_chth': [],
                       'pre_chth': [],
                       'main_xml': [],
                       'out_xml': []}
        :return:
        """
        for macro_type, macros in new_macros.items():
            for macro in macros:
                if macro not in self.macros[macro_type]:
                    self.macros[macro_type].append(macro)

    def _simple_macro_update(self, new_macros, mtype):
        """

        :param new_macros:
        :return:
        """
        for macro in new_macros:
            if macro not in self.macros[mtype]:
                self.macros[mtype].append(macro)

    def sel_prep(self):
        """
        Define the select blob, and then we can make a template out of the dict.
        <option value="normal_yes" selected="true">Yes</option>
        :return:
        """
        sel_blob = []
        for sel in self.blob['options']:
            if self.blob['defaultValue'] == sel['name']:
                sel_blob.append({'value': sel['name'], 'selected': 'true'})
            else:
                sel_blob.append({'value': sel['name'], 'selected': 'false'})

        return sel_blob

    def correct_sci(self, val):
        """
        If the defaultValue contains a scientific notation string, change it to something that looks like a float.
        In GATK json files, these are of the format "1.0E-6"
        Going to convert these instead of using xml validator since the smallest one is reasonable to include as float.
        :return:
        """
        try:
            dividend = float(val.split('E-')[0])
            divisor = float('1' + int(val.split('E-')[1]) * '0')
            return str('%f' % (dividend / divisor))
        except:
            return val

    def _value_correct(self, val):
        """
        Combine different operations necessary to provide a value that is acceptable for the wrapper.
        :return:
        """
        if ',' in val:
            return ''
        elif '[' in val and ']' in val:
            return val.replace('[', '').replace(']', '')
        elif val == 'null':
            return ''
        else:
            return self.correct_sci(val)

    def reblob(self):
        """
        Restructuring the information from the json input file to be used downstream.
        <param name="max_records_in_ram" type="integer" size="10" value="500000" min="1" label="Max Records in RAM" help="When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed." />
        :return:
        """
        xml_out = {'name': self.blob['name'].lstrip('-').replace('-', '_'),
                   'argument': self.blob['name'],
                   'type': self._assign_type(self.type),
                   'label': self.blob['name'].lstrip('-').replace('-', ' ').title(),
                   'optional': self.xml_json_req_map[self.blob['required']],
                   'value': self._value_correct(self.blob['defaultValue']),
                   'format': self.assign_format(),
                   'truevalue': self.blob['name'],
                   'falsevalue': '',
                   'checked': self.blob['defaultValue'],
                   'max': self._assign_min_max(self.blob['maxValue']),
                   'min': self._assign_min_max(self.blob['minValue']),
                   'help': self.blob['summary'],
                   'section': self.section}

        if self.is_req_output or self.is_output:
            xml_out['label'] = self._assign_label(xml_out)
        # if self.in_frmt:
        #     xml_out['format'] = self.in_frmt
        # if self.param_format_map(self.blob['name']):
        #     xml_out['format'] = self.param_format_map(self.blob['name'])
        for key in xml_out.keys():
            #['argument', 'checked', 'falsevalue', 'format', 'help', 'label', 'max', 'min', 'name', 'optional', 'truevalue', 'type', 'value']:
            xml_out[key] = escape(xml_out[key], {'"': '&quot;', "'": '&apos;'})
        return xml_out

    def _assign_label(self, format):
        """
        When we're working with an output, would like to add different information to the label field.
        :return:
        """
        cht_tmpl = self.out_label_tmpl
        return cht_tmpl.substitute(format)

    def _assign_min_max(self, val):
        """
        Fix min and max values.
        :return:
        """
        if val in self.xml_json_num_map:
            return self.xml_json_num_map[val]
        elif val[-2:] == '.0':
            # or int(float(val))
            # or printf-style
            return val[:-2]
        else:
            return val

    def _assign_type(self, type):
        """
        Make adjustments to the 'type' in reblob, for non-standard situations.
        :return:
        """
        if self.is_input:
            return 'data'
        else:
            return type


    def assign_format(self):
        """
        Assign the format where necessary, like if we are dealing with an output param.
        :return:
        """
        if self.is_output:
            if self.pname in self.tool_output_file_type[self.tool_name]:
                return self.tool_output_file_type[self.tool_name][self.pname]
                                       #[self.tool_output_file_type[self.tool_name][self.pname]]
            # if self.pname in self.tool_data[self.tool_name]['output_fmt']:
            #     return self.tool_data[self.tool_name]['output_fmt'][self.pname]
            elif self.pname in self.gen_out_fmt:
                return self.gen_out_fmt[self.pname]
            elif self.pname in self.out_create_params:
                return self.out_create_params[self.pname]
            else:
                return None
        elif self.is_input:
            if self.pname in self.tool_file_type[self.tool_name]:
                return self.file_type_map[self.tool_file_type[self.tool_name][self.pname]]
            # if self.pname in self.tool_data[self.tool_name]['input_fmt']:
            #     return self.tool_data[self.tool_name]['input_fmt'][self.pname]
            elif self.pname in self.gen_in_fmt:
                try:
                    return self.file_type_map[self.gen_in_fmt[self.pname]]
                except:
                    return self.gen_in_fmt[self.pname]
            else:
                return 'txt'
        else:
            # Not sure yet what this will be used for, but I think we need it.
            return ''


class JsonShell(Mappings, XmlTemplates):
    """
    if picard, don't provide macro for annotations, not necessary
    summary == description
    description == help
    name == name
    """
    def __init__(self, args, profile='18.05'):
        # I don't want XmlTemplates to be inherited from multiple places, though right now it is all over the place.  Needs
        # some significant refactoring.
        Mappings.__init__(self)
        XmlTemplates.__init__(self)
        self.profile = profile
        self.args = args
        # Code that appears before the command invocation within the Cheetah section.
        self.pre_chth = [self.chth_tmpl.substitute(macro='set_sections')]
        # Code involved in the actual command invocation of the tool, and parametert
        # We need set_sections for many of the parameters, since we won't know which section these go in until the wrappers
        # are built.  Probably should move this process out of _init_ in to own function.
        self.tool_chth = []
        self.tool_xml = []
        self.xml_opt = []
        self.xml_adv = []
        self.xml_out = []
        self.xml_out_from_work_dir = []
        self.xml_req_out = []
        self.xml_comm = []
        self.xml_dep = []
        self.sel_dict = {}
        self.macros = {'main_chth': {},
                       'pre_chth': {},
                       'main_xml': {},
                       'out_xml': {}
                       }


        with open(args.json, 'rU') as myfile:
            self.json_file = json.load(myfile)
            self.shell_dict = self.build_shell_dict()
            for entry in self.json_file['arguments']:
                self.my_xml = JsonXml(entry, self.shell_dict['short_name'], args)
                if self.my_xml.pname not in self.ignore_params:
                    # If something is in the macros_xml variable, add it here to my_macros.
                    self._macro_update(self.my_xml.macros)

                    # if self.my_xml.pname == 'intervals':
                    #     print(self.my_xml.macros)
                    #     print(self.macros)

                    if self.my_xml.chth:
                        self.tool_chth.extend(self.my_xml.chth)
                    if self.my_xml.chth_pre:
                        self.pre_chth.extend(self.my_xml.chth_pre)

                    if self.my_xml.is_output and not self.my_xml.is_req_output and self.my_xml.pname in self.gen_out_fmt:
                        self.xml_out.append(self.my_xml.xml_param)
                    elif self.my_xml.is_output and not self.my_xml.is_req_output and self.my_xml.pname in self.out_create_params:
                        self.xml_out_from_work_dir.append(self.my_xml.xml_param)
                    elif self.my_xml.is_req_output and not self.my_xml.macros['out_xml']:
                        self.xml_req_out.append(self.my_xml.xml_param)
                    elif self.my_xml.xml_param:
                        if self.my_xml.section == 'advanced':
                            self.xml_adv.append(self.my_xml.xml_param)
                        elif self.my_xml.section == 'optional':
                            self.xml_opt.append(self.my_xml.xml_param)
                        elif self.my_xml.section == 'common':
                            self.xml_comm.append(self.my_xml.xml_param)
                        elif self.my_xml.section == 'deprecated':
                            self.xml_dep.append(self.my_xml.xml_param)
                        elif not self.my_xml.macros['main_xml']:
                            self.tool_xml.append(self.my_xml.xml_param)

                    if self.my_xml.sel_blob:
                        self.sel_dict[self.my_xml.pname] = self.my_xml.sel_blob

    def _macro_update(self, new_macros):
        """
        Update the macros structure.  Structure should look like:
        self.macros = {'main_chth': {<section>: []},
                       'pre_chth': [],
                       'main_xml': [],
                       'out_xml': []}
        :return:
        """
        section = self.my_xml.section
        for macro_type, macros in new_macros.items():
            if section not in self.macros[macro_type]:
                self.macros[macro_type][section] = []
            for macro in macros:
                if macro not in self.macros[macro_type][section]:
                    self.macros[macro_type][section].append(macro)

    def build_shell_dict(self):
        """
        This will house all values the templates need.
        :return:
        """
        shell_dict = {'id': self.json_file['name'].lower().split(' ')[0],
                      'name': Template('GATK4 AUTO $name').substitute(self.json_file),
                      'short_name': self.json_file['name'].split(' ')[0],
                      'profile': self.profile,
                      'description': self.json_file['summary'].rstrip(' '),
                      'summary': pypandoc.convert_text(self.json_file['description'], 'rst', format='html')}
        return shell_dict

    def inputs_create(self):
        """
        Arrange the inputs section.
        :return:
        """
        inputs = []
        for macro in self.my_xml.tool_data[self.shell_dict['short_name']]['pre_params']:
            inputs.append(self.my_xml.xml_tmpl.substitute(macro_name=macro))
        for macro in self.my_xml.tool_data[self.shell_dict['short_name']]['post_params']:
            inputs.append(self.my_xml.xml_tmpl.substitute(macro_name=macro))

        return '\n'.join(inputs)


    def command_create(self):
        """
        Arrange the commands section.
        :return:
        """
        command = []
        command.extend(self.pre_chth)
        command.append(Template('@CMD_BEGIN@ $short_name').substitute(self.shell_dict))
        command.extend(self.tool_chth)
        return '\n'.join(command)


    def xml_chth_expand(self, start_str, in_tmpl, *args, pre=False):
        """
        Based on the json_type, expand the template
        :return:
        """
        temp_str = start_str
        for macros in args:
            for entry in macros:
                if pre:
                    temp_str = in_tmpl.substitute(macro_name=entry)
                    temp_str += start_str
                else:
                    temp_str += in_tmpl.substitute(macro_name=entry)

        return Template(temp_str)


class XmlEtrees(JsonShell):
    """
    Hold all of the etrees we need for the structure.
    """
    def __init__(self, args, id_prefix='gatk4_auto_dev_'):
        """
        Provide templates for the shell of the XML file.
        :return:
        """
        #        etree.write(stdout, xml_declaration=True, encoding='UTF-8')
        JsonShell.__init__(self, args)
        self.args = args
        tool = etree.Element('tool', id=id_prefix + self.shell_dict['id'], name=self.shell_dict['name'],
                             version="@WRAPPER_VERSION@0", profile=self.profile)
        description = etree.SubElement(tool, 'description')
        description.text = '- ' + self.shell_dict['description']
        macros = etree.SubElement(tool, 'macros')
        macros_imp = etree.SubElement(macros, 'import')
        macros_imp.text = 'macros.xml'
        exp_reqs = etree.SubElement(tool, 'expand', macro='requirements')
        exp_vers = etree.SubElement(tool, 'expand', macro='version_cmd')
        command = etree.SubElement(tool, 'command', detect_errors='exit_code')
        command.text = etree.CDATA(self.command_create())

        # INPUT section
        self.inputs = etree.SubElement(tool, 'inputs')

        # Check if there are macros listed in the required section, then place them at the beginning on the xml section.
        if 'required' in self.macros['main_xml']:
            for entry in self.macros['main_xml']['required']:
                etree.SubElement(self.inputs, 'expand', macro=entry)
            self.build_inputs(self.tool_xml, self.inputs, 'param')

        self._section_write('optional', self.xml_opt)
        self._section_write('advanced', self.xml_adv)
        self._section_write('common', self.xml_comm)
        self._section_write('deprecated', self.xml_dep)

        # if self.xml_opt or 'optional' in self.macros['main_xml']:
        #     self.opt_sect = etree.SubElement(self.inputs, 'section', name='optional', title='Optional Parameters', expanded='False')
        #     if 'optional' in self.macros['main_xml']:
        #         for entry in self.macros['main_xml']['optional']:
        #             etree.SubElement(self.opt_sect, 'expand', macro=entry)
        #     self.build_inputs(self.xml_opt, self.opt_sect, 'param')


        # if self.xml_adv or 'advanced' in self.macros['main_xml']:
        #     self.adv_sect = etree.SubElement(self.inputs, 'section', name='advanced', title='Advanced Parameters', expanded='False')
        #     if 'advanced' in self.macros['main_xml']:
        #         for entry in self.macros['main_xml']['advanced']:
        #             etree.SubElement(self.adv_sect, 'expand', macro=entry)
        #     self.build_inputs(self.xml_adv, self.adv_sect, 'param')

        # for entry in self.my_xml.tool_data[self.shell_dict['short_name']]['post_params']:
        #     etree.SubElement(self.inputs, 'expand', macro=entry)

        # # Temporary common section to work on macros.
        # if self.xml_comm or 'common' in self.macros['main_xml']:
        #     self.comm_sect = etree.SubElement(self.inputs, 'section', name='common', title='Common Parameters', expanded='False')
        #     if 'common' in self.macros['main_xml']:
        #         for entry in self.macros['main_xml']['common']:
        #             etree.SubElement(self.comm_sect, 'expand', macro=entry)
        #     self.build_inputs(self.xml_comm, self.comm_sect, 'param')
        #
        #
        # if self.xml_dep:
        #     self.dep_sect = etree.SubElement(self.inputs, 'section', name='deprecated', title='Deprecated Parameters',
        #                                       expanded='False')
        #     self.build_inputs(self.xml_dep, self.dep_sect, 'param')

        # OUTPUT section
        if self.xml_out:
            self.output_sect = etree.SubElement(self.inputs, 'section', name='output_opt',
                                                title='Additional Output Parameters',
                                                expanded='False')
            self.build_inputs_out_sel(self.xml_out, self.output_sect)
            self.build_from_work_dir_out(self.xml_out_from_work_dir, self.output_sect)

        self.outputs = etree.SubElement(tool, 'outputs')
        if 'required' in self.macros['out_xml']:
            for entry in self.macros["out_xml"]['required']:
                etree.SubElement(self.outputs, 'expand', macro=entry)

        self.build_inputs(self.xml_req_out, self.outputs, 'data')
        self.build_inputs(self.xml_out, self.outputs, 'data', True)
        self.build_inputs(self.xml_out_from_work_dir, self.outputs, 'data', True, True)

        tests = etree.SubElement(tool, 'tests')
        help = etree.SubElement(tool, 'help')
        help.text = etree.CDATA(self.shell_dict['summary'])
        citations = etree.SubElement(tool, 'citations')
        exp_cit = etree.SubElement(citations, 'expand', macro='citations')

        self.to_write = etree.tostring(tool, pretty_print=True, encoding="unicode")

    # def _section_write(self, sname, stitle, selname):
    #     """
    #     Write a section, or write a conditional, depending on arg.
    #     :return:
    #     """
    #     if not self.args.old_galaxy:
    #         this_sect = etree.SubElement(self.inputs, 'section', name=sname, title=stitle, expanded='False')
    #         when_yes = None
    #     else:
    #         this_sect = etree.SubElement(self.inputs, 'conditional', name=sname)
    #         this_sect_sel = etree.SubElement(this_sect, 'param', name=selname, type='select',
    #                                                 label=stitle)
    #         opt_yes = etree.SubElement(this_sect_sel, 'option', value='yes')
    #         opt_yes.text = 'yes'
    #         opt_no = etree.SubElement(this_sect_sel, 'option', value='no', selected='true')
    #         opt_no.text = 'no'
    #         when_yes = etree.SubElement(this_sect, 'when', value='yes')
    #     return this_sect, when_yes

    def _section_write(self, section, xml_sect):
        """
        Put section together
        :param section:
        :param xml_sect:
        :return:
        """
        if xml_sect or section in self.macros['main_xml']:
            title = ' '.join([section.title(), 'Parameters'])
            this_sect = etree.SubElement(self.inputs, 'section', name=section, title=title, expanded='False')
            if section in self.macros['main_xml']:
                for entry in self.macros['main_xml'][section]:
                    etree.SubElement(this_sect, 'expand', macro=entry)
            self.build_inputs(xml_sect, this_sect, 'param')

    def _output_section_write(self):
        """
        Write a section, or write a conditional, depending on arg.
        :return:
        """
        if not self.args.old_galaxy:
            self.output_sect = etree.SubElement(self.inputs, 'section', name='output_opt', title='Additional Output Parameters', expanded='False')
        else:
            self.output_sect = etree.SubElement(self.inputs, 'conditional', name='output_opt')
            self.output_sect_sel = etree.SubElement(self.output_sect, 'param', name='output_opt_sel', type='select',
                                                    label='Additional output parameters?')
            self.opt_yes = etree.SubElement(self.output_sect_sel, 'option', value='yes')
            self.opt_yes.text = 'yes'
            self.opt_no = etree.SubElement(self.output_sect_sel, 'option', value='no', selected='true')
            self.opt_no.text = 'no'
            self.when_yes = etree.SubElement(self.output_sect, 'when', value='yes')

    def _build_conditional(self):
        """
        Put the conditional together that acts as a section.
    Example:
    <conditional name="output_opt">
      <param name="output_opt_sel" type="select" label="Additional output parameters?">
        <option value="yes">yes</option>
        <option value="no" selected="true">no</option>
      </param>
      <when value="yes">
        <param argument="--graph_output_sel" checked="false" falsevalue="" help="Write debug assembly graph information to this\
 file" label="Graph Output" name="graph_output_sel" optional="true" truevalue="--graph_output_sel" type="boolean"/>
      </when>
    </conditional>

        :return:
        """
        # HERE
        self.output_sect = etree.SubElement(self.inputs, 'conditional', name='output_opt')
        self.output_sect_sel = etree.SubElement(self.output_sect, 'param', name='output_opt_sel', type='select', label='Additional output parameters?')
        self.opt_yes = etree.SubElement(self.output_sect_sel, 'option', value='yes')
        self.opt_yes.text = 'yes'
        self.opt_no = etree.SubElement(self.output_sect_sel, 'option', value='no', selected='true')
        self.opt_no.text = 'no'
        self.when_yes = etree.SubElement(self.output_sect, 'when', value='yes')


    def build_inputs_out_sel(self, params, parent):
        """
        Build extra options under the additional outputs section to provide booleans for choosing optional outputs.
        <param name="sites-only-vcf-output" argument="--sites-only-vcf-output" type="boolean" truevalue="--sites-only-vcf-output" falsevalue="" optional="true" checked="false" label="Sites Only Vcf Output" help="If true, don&amp;apos;t emit genotype fields when writing vcf file output."/>
        :return:
        """
        for param in params:
            new_name = param['name'] + '_sel'
            new_label = param['name'].replace('_', ' ').title()
            new_arg = '--' + new_name
            this_param = etree.SubElement(parent, 'param', name=new_name, argument=new_arg, type='boolean',
                                          truevalue=new_arg, falsevalue='', optional='true', checked='false',
                                          label=new_label, help=param['help'])

    def build_from_work_dir_out(self, params, parent):
        """
        Work with the output parameters that come in as booleans.  These should be set in the output_opt section
        using the given parameter name, and then set in the output section with a "from_work_dir" value.
        """
        for param in params:
            new_label = param['name'].replace('_', ' ').title()
            new_arg = '--' + param['name'].replace('_', '-')
            this_param = etree.SubElement(parent, 'param', name=param['name'], argument=new_arg, type='boolean',
                                          truevalue=new_arg, falsevalue='', optional='true', checked='false',
                                          label=new_label, help=param['help'])

    def create_output_loc(self):
        """
        Create the output file name for writing, based on input folder.
        :return:
        """
        self.output_name = [self.args.xml_out, 'gatk4_auto_dev_' + self.json_file['name'].lower().split(' ')[0] + '.xml']
        if not self.args.xml_out.endswith('/'):
            return '/'.join(self.output_name)
        else:
            return ''.join(self.output_name)

    def write_me(self):
        """
        Write to file.
        :return:
        """
        handle_out = open(self.create_output_loc(), 'w')
        handle_out.write(self.to_write)
        handle_out.close()

    def build_inputs(self, params, parent, elem, filt_tag=False, from_work_dir=False):
        """

        :return:
        """
        for param in params:
#            if param['optional'] == 'true':
            if 'from_work_dir' not in param and from_work_dir:

                param['from_work_dir'] = 'BLAH'

            this_param = etree.SubElement(parent, elem)
            for k, v in param.items():
                if (k == 'min') and (v == ''):
                    pass
                elif (k == 'max') and (v == ''):
                    pass
                else:
                    this_param.set(k, v)

            if param['name'] in self.sel_dict:
                self.build_sel_opt(this_param, self.sel_dict[param['name']])

            if filt_tag and not from_work_dir:
                # <filter>not gzipped_output</filter>
                filt = etree.SubElement(this_param, 'filter')
                filt.text = 'output_opt[\'' + param['name'] + '_sel\']'

            if filt_tag and from_work_dir:
                filt = etree.SubElement(this_param, 'filter')
                filt.text = 'output_opt[\'' + param['name'] + '\']'

    def build_sel_opt(self, this_param, sel_blob):
        """
        Add select option statements if this is of select type.
        [{'summary': '', 'name': 'ALL'}, {'summary': '', 'name': 'OVERLAPPING_ONLY'}]
        :return:
        """
        for sel in sel_blob:
            this_sel = etree.SubElement(this_param, 'option', selected=sel['selected'], value=sel['value'])
            this_sel.text = sel['value']

class PercentTemplate(Template):
    delimiter = '%'


def main():
    """
    Not including (Picard):
    CREATE_INDEX - Handled by Galaxy internally.
    help - duh
    QUIET - can't think of a good reason to include this, people can just ignore it
    version - Should be included, if anywhere, under the version tag
    :return:
    """
    args = supply_args()
    myshell = XmlEtrees(args)
    myshell.write_me()

if __name__ == "__main__":
    main()
