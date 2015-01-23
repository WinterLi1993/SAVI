#!/usr/bin/env python

import subprocess

def escape_special_char(mystr):
	"""fix string for shell commands by escaping quotes and dollar signs. The idea is, we want to be able to use the echo "cmd" | sh construction"""
	return mystr.replace('"','\\"').replace('$','\\$')

def run_cmd(cmd, bool_verbose, bool_getstdout):
	"""Run system cmd"""

	cmd = "echo \"" + cmd + "\" | sh"

	if (bool_verbose): 
		print(cmd)

	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	proc.wait() 
	# print return code
	# print(proc.returncode) 
	# print stdout stderr tuple
	# proc.communicate()

	(stdout, stderr) =  proc.communicate()

	# if error, print it
	if stderr:
		print("ERROR: " + stderr),

	# return stdout
	if (bool_getstdout): 
		return stdout.rstrip()
	else:
		return "0" # note: this must return a str

def run_qsub(cmd, bool_verbose, bool_getstdout, qsubstr):
	"""Run SGE qsub cmd"""

	cmd = "echo \"" + cmd + "\" | " + qsubstr

	if (bool_verbose): 
		print(cmd)

	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	proc.wait() 

	(stdout, stderr) =  proc.communicate()

	# if error, print it
	if stderr:
		print("ERROR: " + stderr),

	# return job ID, assuming it's in the form:
	# Your job 127748 ("testawesome2") has been submitted
	if (bool_getstdout): 
		return stdout.split()[2]
	else:
		return "0" # note: this must return a str

def whichcmd(cmd, args, wantreturn, wantqsub=0, jobname="myjob", holdstr="0", wantsync=0):
	"""Run cmd as regular cmd or qsub cmd with SGE"""

	if ( args.sge and wantqsub ):
		qsubstr = "qsub "

		if ( holdstr != "0" ):
			qsubstr = qsubstr + "-hold_jid " + holdstr + " " 

		if ( wantsync ):
			qsubstr = qsubstr + "-sync y "

		qsubstr = qsubstr + "-V -N " + jobname + " -e " + args.sgelog + "  -o " + args.sgelog + " -l mem=" + args.sgemem + "G,time=" + args.sgetime + ":: -S /bin/sh -cwd "

		return run_qsub(cmd, args.verbose, wantreturn, qsubstr )
	else: 
		return run_cmd(cmd, args.verbose, wantreturn )
