import unittest
from parse_gatk_json import CheetahPrep

class CheetahPrepTestCase(unittest.TestCase):

    def setUp(self):
        self.pname = "test_arg"
        self.argname = "--test-arg"
        self.mname = "MACRO_TEST"
        self.section = "common"
        self.pre_mname = "PRE_MACRO_TEST"

    def test_req_chth_macro(self):
        self.chth_prep = CheetahPrep(self.pname, self.argname, "required", True, mname=self.mname)
        self.exp_chth = "#include source=$MACRO_TEST#"
        self.assertEqual(self.chth_prep.chth, self.exp_chth)

    def test_req_chth_premacro(self):
        self.chth_prep = CheetahPrep(self.pname, self.argname, "required", True, pre_mname=self.pre_mname)
        self.exp_chth = "#include source=$PRE_MACRO_TEST#"
        self.assertEqual(self.chth_prep.chth, self.exp_chth)

    def test_req_chth(self):
        self.chth_prep = CheetahPrep(self.pname, self.argname, "required", True)
        self.exp_chth = "--test-arg $test_arg"
        self.assertEqual(self.chth_prep.chth, self.exp_chth)

    def test_req_chth_bool(self):
        self.chth_prep = CheetahPrep(self.pname, self.argname, "required", is_req=True, is_bool=True)
        self.exp_chth = "$test_arg"
        self.assertEqual(self.chth_prep.chth, self.exp_chth)

    def test_opt_chth(self):
        self.chth_prep = CheetahPrep(self.pname, self.argname, self.section)
        self.exp_chth = "#if $common.test_arg\n  --test-arg $common.test_arg\n#end if\n"
        self.assertEqual(self.chth_prep.chth, self.exp_chth)

    def test_opt_chth_bool(self):
        self.chth_prep = CheetahPrep(self.pname, self.argname, self.section, is_bool=True)
        self.exp_chth = "#if $common.test_arg\n  $common.test_arg\n#end if\n"
        self.assertEqual(self.chth_prep.chth, self.exp_chth)

    def test_opt_chth_macro(self):
        self.chth_prep = CheetahPrep(self.pname, self.argname, self.section, mname=self.mname)
        self.exp_chth = "#include source=$MACRO_TEST#"
        self.assertEqual(self.chth_prep.chth, self.exp_chth)

    def test_opt_chth_premacro(self):
        self.chth_prep = CheetahPrep(self.pname, self.argname, self.section, pre_mname=self.pre_mname)
        self.exp_chth = "#include source=$PRE_MACRO_TEST#"
        self.assertEqual(self.chth_prep.chth, self.exp_chth)



    # def test_req_chth_macro(self):
    #     self.chth_prep = CheetahPrep(self.pname, "required", True)
    #     self.exp_chth = "#if $%section.%out_sel_name\n%argument $%name\n#end if"
    #     self.assertEqual(self.chth_prep.chth, self.exp_chth)

if __name__ == '__main__':
    unittest.main()



#         self.chth_tmpl = PercentTemplate('#include source=$%macro#')
#         self.req_out_chth = PercentTemplate('%argument $%name')
#         self.vcf_choose = PercentTemplate('#if $%section.%name'
#                                           '\n#if $%section.%name.is_of_type("vcf_bgzip")'
#                                           '\n%argument %name.vcf.gz'
#                                           '\n#else'
#                                           '\n%argument %name.vcf'
#                                           '\n#end if'
#                                           '\n#end if')
#         self.vcf_tabix = PercentTemplate('#if $%section.%name'
#                                          '\n#set datatype = $%section.%name.datatype'
#                                          '\n#if $%section.%name.is_of_type("vcf_bgzip")'
#                                          '\nln -s $%section.%name %name.vcf.gz &&'
#                                          '\ntabix %name.vcf.gz &&'
#                                          '\n#else'
#                                          '\nln -s $%section.%name %name.vcf &&'
#                                          '\n#end if'
#                                          '\n#end if')
#         self.file_chth = PercentTemplate('#if $%section.%out_sel_name\n%argument $%name\n#end if')
# #        self.file_chth_old_gal = PercentTemplate('#if str($output_opt.output_opt_sel) == "yes":\n#if $output_opt.%out_sel_name:\n%argument $%name\n#end if\n#end if')
#         self.ext_arg = '#if $%section.%name\n  %argument $%section.%name\n#end if\n'
# #        self.ext_arg_old_gal = '#if str($%section.%{section}_sel) == "yes":\n#if $%section.%name:\n%argument $%section.%name\n#end if\n#end if'
#         self.reg_arg = '#if $%name\n  %argument $%name\n#end if\n'
