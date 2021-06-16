import json
import os
import time
import logging

def add_to_log(s, level = "info"):
	toprint = time.strftime("%H:%M, %m-%d-%Y:") + s
	if level == "info":
		logging.info(toprint)
	elif level == "debug":
		logging.debug(toprint)
		
def set_access_keys(filepath):
	"""
	Sets the user's access keys to the AWS S3 bucket.

	Assumes this file includes only two uncommented lines with keys:
	export AWS_ACCESS_KEY_ID=<key>
	export AWS_SECRET_ACCESS_KEY=<key>
	"""
	keys = {}
	with open(filepath) as fp:
		for i in fp:
			if len(i.rstrip()) and i[0]!="#":
				cmd = i.rstrip().split(" ")
				assert(cmd[0] == "export")
				cmd.pop(0)
				for i in range(len(cmd)-1,-1,-1):
					if len(cmd[i]) == 0:
						cmd.pop(i)
				pair = "".join(cmd).split("=")
				keys[pair[0]] = pair[1]
	for k in keys:
		os.environ[k] = keys[k]
	return

def load_configs(filename):
    with open(filename) as f: 
        data = f.read()	
    configs = json.loads(data)
    validate_configs(configs)
    return configs

def validate_configs(configs):
	assert(configs["aligner"] == "cellranger")
	return
