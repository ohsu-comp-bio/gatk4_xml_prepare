#!/usr/bin/env python

# Potentially make List[String] types a repeat xml element.  Otherwise, it appears these are
# separated by commas on the command line.
# Annotation options need macros.
# Macros for read filter options?  These have their own json files it seems...
# gVCF not currently set to be different than VCF.

from lxml import etree
from string import Template
from xml.sax.saxutils import escape
import argparse
import json
import pypandoc


VERSION = "0.4.0"

def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('json', help='Input JSON')
    parser.add_argument('--xml_out', help='Output Directory')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


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
    tool_info = JsonTool(args.json)
    myshell = XmlFiller(tool_info)
    XmlWriter(myshell.to_write, args.xml_out, tool_info.tool_name).write_me()


class XmlWriter:
    """

    """
    def __init__(self, to_write, outdir, tool_name, id_prefix='gatk4_auto_'):
        self.to_write = to_write
        self.outdir = outdir
        self.id_prefix = id_prefix
        self.tool_name = tool_name.lower().split(' ')[0]

    def write_me(self):
        handle_out = open(self._create_output_loc(), 'w')
        handle_out.write(self.to_write)
        handle_out.close()

    def _create_output_loc(self):
        """
        Create the output file name for writing, based on input folder.
        :return:
        """
        self.output_name = [self.outdir, self.id_prefix + self.tool_name + '.xml']
        if not self.outdir.endswith('/'):
            return '/'.join(self.output_name)
        else:
            return ''.join(self.output_name)


class XmlEtrees:
    """
    Hold all of the etrees we need for the structure.

    """

    def __init__(self, args):
        """
        Provide templates for the shell of the XML file.
        :return:
        """
        self.args = args
        self.shell_dict = self.args.shell_dict
        self.tool = etree.Element('tool', id=self.args.id_prefix + self.shell_dict['id'],
                                  name=self.shell_dict['name'], version="@WRAPPER_VERSION@0",
                                  profile=self.args.profile)
        description = etree.SubElement(self.tool, 'description')
        description.text = '- ' + self.shell_dict['description']
        macros = etree.SubElement(self.tool, 'macros')
        macros_imp = etree.SubElement(macros, 'import')
        macros_imp.text = 'macros.xml'
        exp_reqs = etree.SubElement(self.tool, 'expand', macro='requirements')
        exp_vers = etree.SubElement(self.tool, 'expand', macro='version_cmd')
        self.command = etree.SubElement(self.tool, 'command', detect_errors='exit_code')
        self.inputs = etree.SubElement(self.tool, 'inputs')
        #self.req_inputs = self._section_write('required')
        self.opt_inputs = self._section_write('optional')
        self.adv_inputs = self._section_write('advanced')
        self.comm_inputs = self._section_write('common')
        self.dep_inputs = self._section_write('deprecated')
        self.output_sect = etree.SubElement(self.inputs, 'section', name='output_opt',
                                            title='Additional Output Parameters', expanded='False')
        self.outputs = etree.SubElement(self.tool, 'outputs')
        tests = etree.SubElement(self.tool, 'tests')
        help = etree.SubElement(self.tool, 'help')
        help.text = etree.CDATA(self.shell_dict['summary'])
        citations = etree.SubElement(self.tool, 'citations')
        exp_cit = etree.SubElement(citations, 'expand', macro='citations')

    def _section_write(self, section):
        """
        Put section together
        :param section:
        :param xml_sect:
        :return:
        """
        title = ' '.join([section.title(), 'Parameters'])
        return etree.SubElement(self.inputs, 'section', name=section, title=title, expanded='False')

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

    def build_sel_opt(self, this_param, sel_blob):
        """
        Add select option statements if this is of select type.
        [{'summary': '', 'name': 'ALL'}, {'summary': '', 'name': 'OVERLAPPING_ONLY'}]
        :return:
        """
        for sel in sel_blob:
            this_sel = etree.SubElement(this_param, 'option', selected=sel['selected'], value=sel['value'])
            this_sel.text = sel['value']


class XmlTemplates(object):
    def __init__(self):
        # Cheetah section
        self.chth_tmpl = PercentTemplate('#include source=$%macro#')
        # Required argument.
        self.req_out_chth = PercentTemplate('%argument $%name')
        # Required boolean argument.
        self.req_out_chth_bool = PercentTemplate('$%name')
        # May want to move the below templates out to the macros section.  These are a little more general
        # since they may apply to any given non-required vcf argument.
        self.vcf_choose = PercentTemplate('#if $%section.%name'
                                          '\n#if $%section.%name.is_of_type("vcf_bgzip")'
                                          '\n%argument %name.vcf.gz'
                                          '\n#else'
                                          '\n%argument %name.vcf'
                                          '\n#end if'
                                          '\n#end if')
        self.vcf_choose_req = PercentTemplate('#if $%name'
                                          '\n#if $%name.is_of_type("vcf_bgzip")'
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
        self.vcf_tabix_req = PercentTemplate('#if $%name'
                                         '\n#set datatype = $%name.datatype'
                                         '\n#if $%name.is_of_type("vcf_bgzip")'
                                         '\nln -s $%name %name.vcf.gz &&'
                                         '\ntabix %name.vcf.gz &&'
                                         '\n#else'
                                         '\nln -s $%name %name.vcf &&'
                                         '\n#end if'
                                         '\n#end if')
        self.file_chth = PercentTemplate('#if $%section.%out_sel_name\n%argument $%name\n#end if')
        self.ext_arg = PercentTemplate('#if $%section.%name\n  %argument $%section.%name\n#end if\n')
        self.ext_arg_bool = PercentTemplate('#if $%section.%name\n  $%section.%name\n#end if\n')
        self.reg_arg = '#if $%name\n  %argument $%name\n#end if\n'

        # XML section. Most of this is handled via etrees, but some cases are more easily handled here.
        self.xml_tmpl = Template('<expand macro="$macro_name" />')
        self.sel_out = Template('<option value="$value" selected="$selected">$value</option>')
        self.out_label_tmpl = Template('$${tool.name} on $${on_string}: $name $format')
        self.out_tmpl = Template('<data format="$format" name="$name" label="$${tool.name} on $${on_string}: $format" />')


class XmlFiller(XmlEtrees, XmlTemplates):
    """

    """
    def __init__(self, tool_args):
        XmlEtrees.__init__(self, tool_args)
        XmlTemplates.__init__(self)
        self.tool_args = tool_args.args
        self.command.text = etree.CDATA(self._command_fill())

        self._fill_macros()
        self._fill_params()
        self._fill_output_params()
        # rm_sect = self._blank_sect_rm()
        self.to_write = etree.tostring(self.tool, pretty_print=True, encoding="unicode")

    def _fill_output_params(self):
        """
        Same as _fill_params, but for the outputs element.
        :return:
        """
        for arg in self.tool_args:
            if not arg.macros['out_xml']:
                if arg.is_output and arg.is_req:
                    self._set_param(arg, self.outputs, 'data')
                elif arg.is_output:
                    this_param = self._set_param(arg, self.outputs, 'data')
                    filt = etree.SubElement(this_param, 'filter')
                    filt.text = 'output_opt[\'' + arg.xml_out['name'] + '_sel\']'

    def _fill_params(self):
        """
        Get all of the parameters placed in the correct sections.
        #     self.build_inputs(self.tool_xml, self.inputs, 'param')
        :return:
        """
        for arg in self.tool_args:
            if not arg.macros['main_xml']:
                if arg.kind == 'required':
                    self._set_param(arg, self.inputs, 'param', first=True)
                elif arg.kind == 'optional':
                    self._set_param(arg, self.opt_inputs, 'param')
                elif arg.kind == 'advanced':
                    self._set_param(arg, self.adv_inputs, 'param')
                elif arg.kind == 'common':
                    self._set_param(arg, self.comm_inputs, 'param')
                elif arg.kind == 'deprecated':
                    self._set_param(arg, self.dep_inputs, 'param')
                elif arg.kind == 'output_opt':
                    self._build_inputs_out_sel(arg.xml_out, self.output_sect)
                else:
                    raise Exception("Unrecognized section type %s." % arg.kind)

    def _build_inputs_out_sel(self, param, parent):
        """
        Build extra options under the additional outputs section to provide booleans for choosing optional outputs.
        <param name="sites-only-vcf-output" argument="--sites-only-vcf-output" type="boolean" truevalue="--sites-only-vcf-output" falsevalue="" optional="true" checked="false" label="Sites Only Vcf Output" help="If true, don&amp;apos;t emit genotype fields when writing vcf file output."/>
        :return:
        """
        new_name = param['name'] + '_sel'
        new_label = param['name'].replace('_', ' ').title()
        new_arg = '--' + new_name
        return etree.SubElement(parent, 'param', name=new_name, argument=new_arg, type='boolean',
                                truevalue=new_arg, falsevalue='', optional='true', checked='false',
                                label=new_label, help=param['help'])

    def _set_param(self, param, parent, elem, first=False):
        """
        :param ToolArgXml
        :return:
        """
        if first:
            this_param = self.inputs.insert(0, etree.SubElement(parent, elem))
        else:
            this_param = etree.SubElement(parent, elem)

        for k, v in param.xml_param.items():
            if (k == 'min') and (v == ''):
                pass
            elif (k == 'max') and (v == ''):
                pass
            else:
                this_param.set(k, v)

        if param.is_sel:
            self._build_sel_opt(this_param, param.sel_blob)

        return this_param

    def _build_sel_opt(self, this_param, sel_blob):
        """
        Add select option statements if this is of select type.
        [{'summary': '', 'name': 'ALL'}, {'summary': '', 'name': 'OVERLAPPING_ONLY'}]
        :return:
        """
        for sel in sel_blob:
            if 'selected' in sel:
                this_sel = etree.SubElement(this_param, 'option', selected=sel['selected'], value=sel['value'])
            else:
                this_sel = etree.SubElement(this_param, 'option', value=sel['value'])
            this_sel.text = sel['value']

    def _fill_macros(self):
        """
        Get the macro tags placed correctly in each section.
        macros = {'main_chth': [],
                  'pre_chth': [],
                  'main_xml': [],
                  'out_xml': []}
        :return:
        """
        for arg in self.tool_args:
            main_macs = arg.macros['main_xml']
            out_macs = arg.macros['out_xml']
            for macro in main_macs:
                if arg.kind == 'required':
                    self.inputs.insert(0, etree.SubElement(self.inputs, 'expand', macro=macro))
                elif arg.kind == 'optional':
                    etree.SubElement(self.opt_inputs, 'expand', macro=macro)
                elif arg.kind == 'advanced':
                    etree.SubElement(self.adv_inputs, 'expand', macro=macro)
                elif arg.kind == 'common':
                    etree.SubElement(self.comm_inputs, 'expand', macro=macro)
                elif arg.kind == 'deprecated':
                    etree.SubElement(self.dep_inputs, 'expand', macro=macro)
                else:
                    raise Exception("Unrecognized section type %s." % arg.kind)
            for macro in out_macs:
                etree.SubElement(self.outputs, 'expand', macro=macro)

    def _blank_sect_rm(self):
        """
        Remove section elements that don't have any children. Provide list of removed elements.
        :return:
        """
        removed = []
        for element in self.inputs:
            if not element.getchildren():
                self.inputs.remove(element)
                removed.append(element)
        return removed

    def _command_fill(self):
        """

        :return:
        """
        command = []
        command.append(self.chth_tmpl.substitute(macro='set_sections'))
        for arg in self.tool_args:
            command.extend(arg.chth_pre)
        command.append(Template('@CMD_BEGIN@ $short_name').substitute(self.shell_dict))
        for arg in self.tool_args:
            command.extend(arg.chth)
        return '\n'.join(command)


class Mappings:
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
                             'GatherType': 'select',
                             'ClippingRepresentation': 'select',
                             'File': 'data',
                             'Set[File]': 'data',
                             'Set[MetricAccumulationLevel]': 'text',
                             'List[File]': 'data',
                             'List[GATKPath]': 'text',
                             'GATKPath': 'text',
                             'List[Double]': 'text',
                             'List[Integer]': 'integer',
                             'List[String]': 'text',
                             'List[Type]': 'text',
                             'List[FeatureInput[Feature]]': 'data',
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


        # For a parameter to be included here, it must have a consistent file format throughout all GATK4 JSONs.
        # Otherwise, this will need to be specified as a macro, and macro will need to be connected to tool.
        self.gen_out_fmt = {
                            'activity_profile_out': 'tabular',
                            'assembly_region_out': 'tabular',
                            'bam_output': 'bam',
                            'graph_output': 'txt',
                            'output_statistics': 'txt'
                            }

        # This is meant to be for those parameters that are booleans, but instruct for an output file to be created.
        # This is a work in progress, will probably ignore these parameters for the time being.
        self.out_create_params = {}

        self.gen_in_fmt = {
            'alleles': 'vcf',
            'allelic_counts': 'tabular',
            'annotated_intervals': 'gatk_interval',
            'BAIT_INTERVALS': 'picard_interval_list',
            'bam_shard': 'data',
            'bam_with_header': 'sam',
            'clip_sequences_file': 'fasta',
            'comp': 'vcf',
            'comparison': 'vcf',
            'concordance': 'vcf',
            'contamination_fraction_per_sample_file': 'tabular',
            'count_panel_of_normals': 'h5',
            'dbsnp': 'vcf',
            'denoised_copy_ratios': 'tabular',
            'discordance': 'vcf',
            'feature_file': 'vcf',
            'gatk_config_file': 'txt',
            'germline_resource': 'vcf',
            'intervals': 'gatk_interval',
            'known_sites': 'vcf',
            'normal_allelic_counts': 'tabular',
            'pedigree': 'tabular',
            'population_callset': 'vcf',
            'panel_of_normals': 'vcf',
            'read_index': 'idx',
            'segments': 'tabular',
            'TARGET_INTERVALS': 'picard_interval_list',
            'variant': 'vcf'
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

        # For case where there is no required output parameter (output goes to stdout).
        self.stdout_tools = ['CountBases', 'CountReads', 'FlagStat']
        self.stdout_chth = ['stdout_to_output']
        self.stdout_macros = ['stdout_to_output_params']

        # Define which fields are associated with which xml data types.
        self.param_tmpls = {'integer': ['name', 'argument', 'type', 'optional', 'value', 'min', 'max', 'label', 'help'],
                            'float': ['name', 'argument', 'type', 'optional', 'value', 'min', 'max', 'label', 'help'],
                            'text': ['name', 'argument', 'type', 'optional', 'value', 'label', 'help'],
                            'data': ['name', 'argument', 'type', 'optional', 'format', 'label', 'help'],
                            'select': ['name', 'argument', 'type', 'optional', 'label', 'help'],
                            'boolean': ['name', 'argument', 'type', 'truevalue', 'falsevalue', 'optional', 'checked', 'label', 'help'],
                            'output': ['format', 'name', 'label']}

        # Many parameters will take multiple Galaxy datatypes.  Providing a mapping will give flexibility to add or
        # subtract types when wanted.
        self.file_type_map = {'vcf': 'vcf,vcf_bgzip',
                              'idx': 'tabix,bai',
                              'gatk_interval': 'gatk_interval,bed,vcf',
                              'sam': 'sam,bam',
                              'picard_interval_list': 'picard_interval_list,tabular',
                              'tabular': 'tabular'}

        # Relate tool name to a common argument name, such as input (which can be connected to different file types),
        # and a file type.
        self.tool_file_type = {'AnnotatePairOrientation': {'input': 'sam'},
                               'BaseRecalibrator': {'input': 'sam'},
                               'BwaMemIndexImageCreator': {'input': 'fasta'},
                               'ClipReads': {'input': 'sam'},
                               'CombineGVCFs': {'variant': 'vcf',
                                                'input': 'sam'},
                               'CountBases': {'input': 'sam'},
                               'CountReads': {'input': 'sam'},
                               'FixCallSetSampleOrdering': {'input': 'sam'},
                               'FixMisencodedBaseQualityReads': {'input': 'sam'},
                               'FlagStat': {'input': 'sam'},
                               'GatherVcfsCloud': {'input': 'vcf'},
                               'GetSampleName': {'input': 'sam'},
                               'HaplotypeCaller': {'input': 'sam'},
                               'LeftAlignIndels': {'input': 'sam'},
                               'Mutect2': {'input': 'sam'},
                               'PrintReads': {'input': 'sam'},
                               'SortSam': {'INPUT': 'sam'},
                               'SplitReads': {'input': 'sam'},
                               'IntervalListToBed': {'INPUT': 'picard_interval_list'}}

        self.tool_output_file_type = {'AnnotatePairOrientation': {'output': 'vcf'},
                                      'BaseRecalibrator': {'output': 'gatk_recal'},
                                      'BwaMemIndexImageCreator': {'output': 'data'},
                                      'ClipReads': {'output': 'sam'},
                                      'CombineGVCFs': {'output': 'vcf'},
                                      'ConvertHeaderlessHadoopBamShardToBam': {'output': 'sam'},
                                      'FixCallSetSampleOrdering': {'output': 'vcf'},
                                      'FixMisencodedBaseQualityReads': {'output': 'sam'},
                                      'GatherVcfsCloud': {'output': 'vcf'},
                                      'GetSampleName': {'output': 'txt'},
                                      'HaplotypeCaller': {'output': 'vcf'},
                                      'IndexFeatureFile': {'output': 'data'},
                                      'LeftAlignIndels': {'OUTPUT': 'sam'},
                                      'Mutect2': {'output': 'vcf'},
                                      'PrintReads': {'output': 'sam'},
                                      'SortSam': {'OUTPUT': 'sam'},
                                      'IntervalListToBed': {'OUTPUT': 'bed'}}

        # Relate the name of the argument to the selection of macros that goes with it.
        self.macro_to_param = {'input': {'sam': {'pre_chth': ['bam_index_pre_chth'],
                                                 'main_chth': ['gatk_bam_input'],
                                                 'main_xml': ['gatk_bam_req_params']},
                                         'fasta': {'main_chth': ['ref_opts_input'],
                                                   'main_xml': ['ref_sel']},
                                         'vcf': {'pre_chth': [],
                                                 'main_chth': [],
                                                 'main_xml': []}
                                         },
                               'variant': {'vcf': {'pre_chth': ['gatk_tabix_multi'],
                                                   'main_chth': ['gatk_input_multi'],
                                                   'main_xml': ['vcf_input_params_multi']}
                                           },
                               'INPUT': {'sam': {'pre_chth': ['bam_index_pre_chth'],
                                                 'main_chth': ['picard_bam_input'],
                                                 'main_xml': ['gatk_bam_req_params']}
                                         },
                               'output': {'vcf': {'main_chth': ['vcf_output_opts'],
                                                  'main_xml': ['gzip_vcf_params'],
                                                  'out_xml': ['gzip_vcf_output_params']}
                                          }
                               }


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
    def __init__(self, pname, argname, section, is_req=False, is_input_vcf=False, mname=None, pre_mname=None,
                 is_bool=False, out_sel_name=None):
        XmlTemplates.__init__(self)
        # Output is required, and there is a macro in my_xml.macros_tmpl.

        self.pname = pname
        self.argname = argname
        self.section = section
        self.mname = mname
        self.pre_mname = pre_mname
        self.is_req = is_req
        self.is_input_vcf = is_input_vcf
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
        if self.is_input_vcf and not self.is_req and not self.pre_mname:
            return self.vcf_choose.substitute(section=self.section, name=self.pname, argument=self.argname)

        if self.is_input_vcf and self.is_req and not self.pre_mname:
            return self.vcf_choose_req.substitute(section=self.section, name=self.pname, argument=self.argname)

        if self.is_input_vcf and not self.is_req and self.pre_mname:
            return self.vcf_tabix.substitute(section=self.section, name=self.pname)

        if self.is_input_vcf and self.is_req and self.pre_mname:
            return self.vcf_tabix_req.substitute(name=self.pname)

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

        raise Exception("NO CHEETAH TEMPLATE DECLARED.")


class JsonTool:
    """
    Contains tool-level information, from first-level of tool json.
    In general, fields are:
    summary
    arguments
    description
    name
    group
    beta
    experimental

    org_broadinstitute_hellbender_tools_PrintReads.json contains additional field (why):
    walkerType

    Will expose what we need and ignore the rest.  arguments will be broken up in to separate ToolArg classes.
    """
    def __init__(self, filename, profile='18.01', id_prefix='gatk4_auto_'):
        self.filename = open(filename, 'r')
        self.profile = profile
        self.id_prefix = id_prefix
        self.tool_desc = json.load(self.filename)

        self.summary = self.tool_desc['summary']
        self.description = self.tool_desc['description']
        self.tool_name = self.tool_desc['name']
        self.shell_dict = self._build_shell_dict()
        self._params_stdout()
        self.ignore = self._ignore_list()
        self.args = self._prep_args()

    def _ignore_list(self):
        """
        Provide a list of parameter names we don't want to include in the wrapper.
        :return:
        """
        # For parameters we universally do not want to see in the UI.
        # interval_padding is taken care of any time we see an intervals parameter.
        return ('help', 'version', 'showHidden', 'interval_padding', 'interval_exclusion_padding', 'TMP_DIR',
                'create_output_bam_index', 'create_output_bam_md5', 'read_index', 'create_output_variant_index',
                'create_output_variant_md5')

    def _build_shell_dict(self):
        """
        This will house all values the templates need.
        :return:
        """
        shell_dict = {'id': self.tool_desc['name'].lower().split(' ')[0],
                      'name': Template('GATK4 AUTO $name').substitute(self.tool_desc),
                      'short_name': self.tool_desc['name'].split(' ')[0],
                      'profile': self.profile,
                      'description': self.tool_desc['summary'].rstrip(' '),
                      'summary': pypandoc.convert_text(self.tool_desc['description'], 'rst', format='html')}
        return shell_dict

    def _prep_args(self):
        """
        Take the arguments portion of tool_desc and make in to ToolArg classes.
        :return:
        """
        args = []
        for arg in self.tool_desc['arguments']:
            name = arg['name'].lstrip('-').replace('-', '_')
            if name not in self.ignore:
                args.append(ToolArgXml(arg, self.tool_name))
        return args

    def _params_stdout(self):
        """
        Print variable contents to stdout.
        :return:
        """
        print("Tool Name: {0}".format(self.tool_name))


class ToolArgBasic:
    """
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

    Set what we can here, without pulling from mappings or templates.  We can deal with rest of mappings with
    superclass.
    """
    def __init__(self, arg):
        self.arg = arg

        self.kind = self.arg['kind']
        self.type = self.arg['type']
        self.arg_name = self.arg['name']
        self.pname = self.arg['name'].lstrip('-').replace('-', '_')

        self.is_bool = self._bool_set()
        self.is_sel = self._sel_set()

        self.sel_blob = self._sel_prep()


    def _bool_set(self):
        """
        Provide convenience variable to denote boolean.
        :return:
        """
        if self.type == 'boolean':
            return True
        return False

    def _sel_set(self):
        """

        :return:
        """
        if self.arg['options']:
            return True
        return False

    def _sel_prep(self):
        """
        Define the select blob, and then we can make a template out of the dict.
        <option value="normal_yes" selected="true">Yes</option>
        :return:
        """
        sel_blob = []
        for sel in self.arg['options']:
            if self.arg['defaultValue'] == sel['name']:
                sel_blob.append({'value': sel['name'], 'selected': 'true'})
            else:
                sel_blob.append({'value': sel['name']})
        return sel_blob


    def _params_stdout(self):
        """
        Print variable contents to stdout.
        :return:
        """
        print("Parameter Name: {0}".format(self.pname))
        print("JSON Type: {0}".format(self.type))
        print("JSON Section: {0}".format(self.kind))
        print("Is Boolean: {0}".format(self.is_bool))
        print("Is Select: {0}".format(self.is_sel))
        print("Select Blob: {0}".format(self.sel_blob))


class ToolArg(ToolArgBasic, Mappings):
    """

    """

    def __init__(self, arg, tool_name):
        """

        """
        # Grab mappings, so that we can provide tool values.
        ToolArgBasic.__init__(self, arg)
        Mappings.__init__(self)

        self.tool_name = tool_name
        self.is_output = self._output_set()
        self.xml_type = self._type_set()
        self.is_input = self._input_set()
        if self.is_input:
            self.xml_type = 'data'
        self.is_input_vcf = self._input_vcf_set()
        self.is_req = self._req_set()
        self.ftype = self._ftype_set()
        self.macros = self._get_macros()
        self._sel_name_set()

    def _sel_name_set(self):
        """

        :return:
        """
        if self.is_output and not self.is_req and self.pname in self.gen_out_fmt:
            self.kind = 'output_opt'
            self.out_sel_name = self.pname + '_sel'
            self.out_sel_arg = self.arg['name'] + '_sel'
        elif self.is_output and not self.is_req and self.pname in self.out_create_params:
            self.kind = 'output_opt'
            self.out_sel_name = self.pname
            self.out_sel_arg = self.arg['name']
        else:
            self.out_sel_name = None
            self.out_sel_arg = None


    def _ftype_set(self):
        """
        Set the file type, for data parameters.
        self.tool_file_type = {'AnnotatePairOrientation': {'input': 'sam'},
                               'BaseRecalibrator': {'input': 'sam'},
                               'BwaMemIndexImageCreator': {'input': 'fasta'},
        self.tool_output_file_type = SAME
        :return:
        """
        if self.pname in self.gen_in_fmt:
            return self.gen_in_fmt[self.pname]
        if self.pname in self.gen_out_fmt:
            return self.gen_out_fmt[self.pname]
        if self.tool_name in self.tool_file_type:
            if self.pname in self.tool_file_type[self.tool_name]:
                return self.tool_file_type[self.tool_name][self.pname]
        if self.tool_name in self.tool_output_file_type:
            if self.pname in self.tool_output_file_type[self.tool_name]:
                return self.tool_output_file_type[self.tool_name][self.pname]
        return 'txt'

    def _get_macros(self):
        """
        Update the macros structure.  Structure should look like:
        self.macros = {'main_chth': [],
                       'pre_chth': [],
                       'main_xml': [],
                       'out_xml': []}
        :return:
        These go in main_xml.
        self.param_to_macro_xml = {'exclude_intervals': ['gatk_excl_ints'],
                                   'intervals': ['gatk_ints']}
        These go in main_chth.
        self.param_to_macro_tmpl = {'exclude_intervals': ['gatk_excl_ints_chth'],
                                    'intervals': ['gatk_ints_chth']}
        These go in pre_chth.
        self.param_to_macro_pretmpl = {'exclude_intervals': ['pre_gatk_excl_ints_chth'],
                                       'intervals': ['pre_gatk_ints_chth']}
        These get placed accordingly.
        self.macro_to_param = {'input': {'sam': {'pre_chth': ['bam_index_pre_chth'],
                                                 'main_chth': ['gatk_bam_input'],
                                                 'main_xml': ['gatk_bam_req_params']},
                                         'fasta': {'main_chth': ['ref_opts_input'],
                                                   'main_xml': ['ref_sel']},
                                         'vcf': {'pre_chth': [],
                                                 'main_chth': [],
                                                 'main_xml': []}
                                         }}
        These will go in to main_chth and out_xml, respectively, and depend on tool name for placement.
        self.stdout_tools = ['CountBases', 'CountReads', 'FlagStat']
        self.stdout_chth = ['stdout_to_output']
        self.stdout_macros = ['stdout_to_output_params']
        """
        macros = {'main_chth': [],
                  'pre_chth': [],
                  'main_xml': [],
                  'out_xml': []}
        if self.pname in self.param_to_macro_xml:
            macros['main_xml'].extend(self.param_to_macro_xml[self.pname])
        if self.pname in self.param_to_macro_tmpl:
            macros['main_chth'].extend(self.param_to_macro_tmpl[self.pname])
        if self.pname in self.param_to_macro_pretmpl:
            macros['pre_chth'].extend(self.param_to_macro_pretmpl[self.pname])
        if self.pname in self.macro_to_param:
            if self.ftype in self.macro_to_param[self.pname]:
                for zone, macro in self.macro_to_param[self.pname][self.ftype].items():
                    macros[zone].extend(macro)
        if self.tool_name in self.stdout_tools:
            macros['main_chth'].extend(self.stdout_chth)
            macros['main_xml'].extend(self.stdout_macros)
        return macros

    def _type_set(self):
        """
        Look at raw type value from JSON, then use the mapping to convert to XML type.
        :return:
        """
        if self.type in self.xml_json_type_map:
            return self.xml_json_type_map[self.type]
        else:
            raise Exception('Argument type %s not recognized.' % self.type)

    def _output_set(self):
        """
        Decide whether we are denoting this arg as output.
        :return:
        """
        if self.pname in self.gen_out_fmt:
            return True
        if self.tool_name in self.tool_output_file_type:
            if self.pname in self.tool_output_file_type[self.tool_name]:
                return True
        if self.pname in self.out_create_params:
            return True
        return False

    def _input_set(self):
        """
        Decide whether we are denoting this arg as input.
        :return:
        """
        if self.pname in self.gen_in_fmt:
            return True
        if self.tool_name in self.tool_file_type:
            if self.pname in self.tool_file_type[self.tool_name]:
                return True
        if self.xml_type == 'data' and not self.is_output:
            return True
        return False

    def _input_vcf_set(self):
        """
        Decide whether we are denoting this arg as input vcf.
        :return:
        """
        if self.pname in self.gen_in_fmt:
            if self.is_input and self.gen_in_fmt[self.pname] == 'vcf':
                return True
        return False

    def _req_set(self):
        """
        Decide whether we are denoting this arg as being required.
        If this arg is in Mappings.macro_to_param, it will not be denoted as required, though in the macro
        it may actually be required.  This is mostly relevant for output args, since input arg macros still will
        need to be placed in the appropriate section.
        :return:
        """
        if self.kind == 'required':
            return True
        if self.pname in self.macro_to_param:
            self.kind = 'required'
            return True
        return False

    def _params_stdout(self):
        """
        Print variable contents to stdout.
        :return:
        """
        super()._params_stdout()
        print("Is Output: {0}".format(self.is_output))
        print("Is Input: {0}".format(self.is_input))
        print("Is Input VCF: {0}".format(self.is_input_vcf))
        print("Is Required: {0}".format(self.is_req))
        print("XML Type: {0}".format(self.xml_type))
        print("File Type: {0}".format(self.ftype))
        print("Macros: {0}".format(self.macros))


class ToolArgXml(ToolArg):
    """
    Here we provide:
    xml_out dict (provided for tool param templating)
    """
    def __init__(self, arg, tool_name):
        """

        """
        ToolArg.__init__(self, arg, tool_name)
        self.xml_out = self.reblob()

        # for mtype, macros in self.macros.items():
        #     for macro in macros:
        #         print(macro)
        #         print(CheetahPrep(self.pname, self.xml_out['argument'], self.kind, self.is_req, self.is_input_vcf,
        #                     mname=macro, is_bool=self.is_bool, out_sel_name=None).chth)

        # Create Cheetah strings.
        self.chth = []
        if self.macros['main_chth']:
            for macro in self.macros['main_chth']:
                self.chth.append(CheetahPrep(self.pname, self.xml_out['argument'], self.kind, self.is_req, self.is_input_vcf,
                                             mname=macro, is_bool=self.is_bool, out_sel_name=self.out_sel_name).chth)
        elif not self.is_input_vcf:
            self.chth.append(CheetahPrep(self.pname, self.xml_out['argument'], self.kind, self.is_req, self.is_input_vcf,
                                         mname=None, is_bool=self.is_bool, out_sel_name=self.out_sel_name).chth)

        self.chth_pre = []
        if self.macros['pre_chth']:
            for macro in self.macros['pre_chth']:
                self.chth_pre.append(CheetahPrep(self.pname, self.xml_out['argument'], self.kind, self.is_req,
                                                 self.is_input_vcf, pre_mname=macro, is_bool=self.is_bool).chth)

        if self.is_input_vcf:
            self.chth.append(CheetahPrep(self.pname, self.xml_out['argument'], self.kind, self.is_req,
                                             self.is_input_vcf).chth)
            self.chth_pre.append(CheetahPrep(self.pname, self.xml_out['argument'], self.kind, self.is_req,
                                             self.is_input_vcf, pre_mname='blah').chth)

        self.xml_param = self._set_xml_param()

    def _set_xml_param(self):
        """
        Get the fields we want per XML type.
        :return:
        """
        if self.is_output:
            xml_param = {label: self.xml_out[label] for label in self.param_tmpls['output']}
        else:
            xml_param = {label: self.xml_out[label] for label in self.param_tmpls[self.xml_out['type']]}
        return xml_param

    def _assign_optional(self, val):
        """

        :return:
        """
        if val == 'no':
            return 'true'
        elif val == 'yes':
            return 'false'
        else:
            raise Exception('JSON required value %s not recognized.' % val)

    def _assign_min_max(self, val):
        """
        self.xml_json_num_map = {'Infinity': '', '-Infinity': '', 'NA': ''}
        Fix min and max values.
        :return:
        """
        bad_nums = ['Infinity', '-Infinity', 'NA']
        if val in bad_nums:
            return ''
        if val[-2:] == '.0':
            return val[:-2]
        return val

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

    def assign_format(self):
        """
        Assign the format where necessary, like if we are dealing with an output param.
        :return:
        """
        if self.is_output:
            if self.tool_name in self.tool_output_file_type:
                if self.pname in self.tool_output_file_type[self.tool_name]:
                    return self.tool_output_file_type[self.tool_name][self.pname]
                elif self.pname in self.gen_out_fmt:
                    return self.gen_out_fmt[self.pname]
                elif self.pname in self.out_create_params:
                    return self.out_create_params[self.pname]
                else:
                    return ''

        if self.is_input:
            if self.tool_name in self.tool_file_type:
                if self.pname in self.tool_file_type[self.tool_name]:
                    try:
                        return self.file_type_map[self.tool_file_type[self.tool_name][self.pname]]
                    except:
                        return self.tool_file_type[self.tool_name][self.pname]
            elif self.pname in self.gen_in_fmt:
                try:
                    return self.file_type_map[self.gen_in_fmt[self.pname]]
                except:
                    return self.gen_in_fmt[self.pname]
            else:
                return 'txt'

    # def _assign_label(self, format):
    #     """
    #     When we're working with an output, would like to add different information to the label field.
    #     :return:
    #     """
    #     cht_tmpl = self.out_label_tmpl
    #     return cht_tmpl.substitute(format)

    def reblob(self):
        """
        Restructuring the information from the json input file to be used downstream.
        <param name="max_records_in_ram" type="integer" size="10" value="500000" min="1" label="Max Records in RAM" help="When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed." />
        :return:
        """
        xml_out = {'name': self.pname,
                   'argument': self.arg_name,
                   'type': self.xml_type,
                   'label': self.arg_name.lstrip('-').replace('-', ' ').title(),
                   'optional': self._assign_optional(self.arg['required']),
                   'value': self._value_correct(self.arg['defaultValue']),
                   'format': self.ftype,
                   'truevalue': self.arg_name,
                   'falsevalue': '',
                   'checked': self.arg['defaultValue'],
                   'max': self._assign_min_max(self.arg['maxValue']),
                   'min': self._assign_min_max(self.arg['minValue']),
                   'help': self.arg['summary'],
                   'section': self.kind}

        # if self.is_output:
        #     xml_out['label'] = self._assign_label(xml_out)
        for key in xml_out:
            # ['argument', 'checked', 'falsevalue', 'format', 'help', 'label', 'max', 'min', 'name', 'optional', 'truevalue', 'type', 'value']:
            if xml_out[key]:
                xml_out[key] = escape(xml_out[key], {'"': '&quot;', "'": '&apos;'})
            else:
                xml_out[key] = ''

        return xml_out


class PercentTemplate(Template):
    delimiter = '%'


if __name__ == "__main__":
    main()
