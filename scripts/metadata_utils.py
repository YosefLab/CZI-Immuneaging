import os
import time
import pandas as pd
import warnings

def read_google_sheet(spreadsheet_id, sheet_name):
	"""
	Downloads and reads a specific sheet from a Google Spreadsheet.
	Parameters
	----------
	spreadsheet_id
		id of a publicly accessible Google Spreadsheet.
	sheet_name
		Specific sheet to return.
	"""
	try:
		import gdown
	except ModuleNotFoundError as e:
		raise ModuleNotFoundError(
			"gdown is not installed. Please install gdown via: pip install gdown"
		)
	
	# url needs to be in a specific format
	with warnings.catch_warnings(record=True) as w:
		warnings.filterwarnings("error")
		url = "https://docs.google.com/spreadsheets/d/{0}/gviz/tq?tqx=out:csv&sheet={1}".format(spreadsheet_id, sheet_name)
		# use a tmp file name since calling gdown.download with output = None still saves a file; remove file at the end
		output_fn = str(time.time()) + ".csv"
		output_fn = gdown.download(url, output = output_fn, quiet=False)
		if len(w) == 1:
			warnings.showwarning(
				msg.message, msg.category, msg.filename, msg.lineno, msg.line
			)
	
	# read spreadsheet
	data = pd.read_csv(output_fn)
	os.remove(output_fn)
	return data
