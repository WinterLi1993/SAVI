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
		mysteps = ['step1', 'step2', 'step3', 'step4', 'step5']

		for i in mysteps: 
			stepobject = getattr(savi_dev, i, None)  
			self.assertTrue(callable(getattr(stepobject, "run", None)))

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

if __name__ == "__main__":
	unittest.main()
