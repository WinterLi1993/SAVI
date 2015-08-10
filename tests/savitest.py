#!/usr/bin/env python

import unittest
import sys

sys.path.append("..") 

from SAVI import savi_dev 

"""
	savitest.py
	~~~~~~
	Unit tests for SAVI 
"""

class saviTests(unittest.TestCase):

	def testForMethodsExistence(self):
		"""Check if run methods exist in step objects"""

		# the steps of SAVI
		mysteps = ['Step1', 'Step2', 'Step3', 'Step4', 'Step5']

		for i in mysteps: 
			stepclass = getattr(savi_dev, i, None)  
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
			self.assertEqual(output, savi_dev.generate_compsamp(input))

	def testGetUniq(self):
		"""Test the get_uniq_samples function - should return list of uniq of samples from the sample comparison string"""

		# expected input output
		inout = ( 
				('2:1,3:1,3:2', ['1', '2', '3']),
			)

		for input, output in inout: 
			self.assertEqual(output, savi_dev.get_uniq_samples(input))

	def testPriorStr(self):
		"""Test the generate_priorstr function - should return comma delim string of index:path"""

		self.assertEqual(savi_dev.generate_priorstr('2:1,3:1,3:2','/test/path'), '1:/test/path,2:/test/path,3:/test/path')

	def testDefaultPrior(self):
		"""Test that diploid prior is used by default"""

		# note that this messes up the arguments for the test script
		(args, parser) = savi_dev.get_arg()
		# print(myargs.priorstring)
		assert 'prior_diploid01' in args.priorstring

if __name__ == "__main__":

	# unittest.main()

	# hardwire verbose == on option (http://stackoverflow.com/questions/13034207/unittest-increase-modules-verbosity-when-tested)
	suite = unittest.TestLoader().loadTestsFromTestCase(saviTests)
	unittest.TextTestRunner(verbosity=2).run(suite)
