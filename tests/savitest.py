#!/usr/bin/env python

import unittest
import sys

sys.path.append("..") 

from SAVI import savi 
import os 
import filecmp

"""
	savitest.py
	~~~~~~
	Unit tests for SAVI 
"""

class saviSimpleFuntionsTests(unittest.TestCase):
	"""Test simple helper functions in main wrapper"""

	def testForMethodsExistence(self):
		"""Check if run methods exist in step objects"""

		# the steps of SAVI
		mysteps = ['Step1', 'Step2', 'Step3', 'Step4', 'Step5']

		for i in mysteps: 
			stepclass = getattr(savi, i, None)  
			self.assertTrue(callable(getattr(stepclass, "run", None)))

	def testCompSamp(self):
		"""Test the sample comparison string: everything should be compared to sample 1 (assumed to be "normal") by default"""

		# expected input output
		inout = ( 
				(1, '1'), 
				(2, '2:1'), 
				(3, '2:1,3:1'), 
				(4, '2:1,3:1,4:1')
			)

		for input, output in inout: 
			self.assertEqual(output, savi.generate_compsamp(input))

	def testGetUniq(self):
		"""Test the get_uniq_samples function - should return list of uniq of samples from the sample comparison string"""

		# expected input output
		inout = ( 
				('2:1,3:1,3:2', ['1', '2', '3']),
			)

		for input, output in inout: 
			self.assertEqual(output, savi.get_uniq_samples(input))

	def testPriorStr(self):
		"""Test the generate_priorstr function - should return comma delim string of index:path"""

		self.assertEqual(savi.generate_priorstr('2:1,3:1,3:2','/test/path'), '1:/test/path,2:/test/path,3:/test/path')

class saviArgsTests(unittest.TestCase):
	"""Tests related to input arguments"""

	def testDefaultPrior(self):
		"""Test that diploid prior is used by default"""

		# note that this messes up the arguments for the test script
		(args, parser) = savi.get_arg()
		# print(myargs.priorstring)
		assert 'prior_diploid01' in args.priorstring

class saviReportsTests(unittest.TestCase):
	"""Tests related to report filtering"""

	def setUp(self):
		"""Set Up"""

		# set attributes
		# directory where this script resides 			
		software = os.path.dirname(os.path.realpath(__file__))
		self.software = software
		# reports 
		self.report_all = self.software + "/test_input_output/report.all.vcf"
		self.report_coding = self.software + "/test_input_output/report.coding.vcf"
		# tmp file
		self.report_out = self.software + "/test_input_output/tmp.vcf"

		# get arguments dictionary
		(args, parser) = savi.get_arg()

		# artificially update the args dict to find finalsavi.vcf, so it doesn't complain
		vars(args)['reportdir'] = self.software + "/test_input_output" 

		# instantiate a step 5 obj 
		self.step5 = savi.Step5(args, "5")

	def tearDown(self):
		"""Clean Up"""

		# clean up
		if os.path.isfile(self.report_out):
			os.remove(self.report_out)

	def testCodingFilter(self):
		"""Test that Coding Filter function works"""

		# run the function
		self.step5.filterCoding(self.report_all, self.report_out)

		# test that output matches self.report_coding
		self.assertTrue(filecmp.cmp(self.report_out, self.report_coding))

	def testNFitler(self):
		"""Test that remove Ns function works"""

		# run the function
		self.step5.filterNs(self.software + "/test_input_output/report.all.add_Ns.vcf", self.report_out)

		# test that output has Ns removed
		self.assertTrue(filecmp.cmp(self.report_out, self.software + "/test_input_output/report.all.remove_Ns.vcf"))

	def testPDStr(self):
		"""Test the PD filter string"""

		# expected input output
		inout = ( 
				('1', ('', '')), 
				('1,2', ('', '')), 
				('2:1', ('PD21_L > 0', 'PD21_U < 0')), 
				('2:1,3:1,5', ('PD21_L > 0 | PD31_L > 0', 'PD21_U < 0 | PD31_U < 0')), 
				('2:1,3:1,3:2', ('PD21_L > 0 | PD31_L > 0 | PD32_L > 0', 'PD21_U < 0 | PD31_U < 0 | PD32_U < 0'))
			)

 		for input, output in inout: 
 			self.assertEqual(output, self.step5.generate_PD_str(input))

if __name__ == "__main__":

	# unittest.main()

	# hardwire verbose == on option (http://stackoverflow.com/questions/13034207/unittest-increase-modules-verbosity-when-tested)
	for i in [saviSimpleFuntionsTests, saviArgsTests, saviReportsTests]:
		suite = unittest.TestLoader().loadTestsFromTestCase(i)
		unittest.TextTestRunner(verbosity=2).run(suite)
