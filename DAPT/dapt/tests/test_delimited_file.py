"""
    Test if delimited_file.py is working correctly
"""

import dapt
import os, csv
from collections import OrderedDict 

def create_test_file():
	# Reset csv.  This is just used update the csv and reset it so interesting things happen.
	with open('test.csv', 'w') as f:
		writer = csv.DictWriter(f, fieldnames=['id', 'startTime', 'endTime', 'status', 'a', 'b', 'c'])
		writer.writeheader()
		writer.writerow({'id':'t1', 'startTime':'2019-09-06 17:23', 'endTime':'2019-09-06 17:36', 'status':'finished', 'a':'2', 'b':'4', 'c':'6'})
		writer.writerow({'id':'t2', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'10', 'c':''})
		writer.writerow({'id':'t3', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'-10', 'c':''})

# Test that we can read a delimited file
def test_DF_read():
	create_test_file()

	db = dapt.delimited_file.Delimited_file('test.csv', ',')
	
	actual = db.get_table()

	expected = []
	expected.append(OrderedDict({'id':'t1', 'startTime':'2019-09-06 17:23', 'endTime':'2019-09-06 17:36', 'status':'finished', 'a':'2', 'b':'4', 'c':'6'}))
	expected.append({'id':'t2', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'10', 'c':''})
	expected.append({'id':'t3', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'-10', 'c':''})

	assert actual == expected, "Cannot read the delimited file.  Nothing else should work."

# Test if the keys from a delimited file can be retrieved
def test_DF_get_keys():
	create_test_file()

	db = dapt.delimited_file.Delimited_file('test.csv', ',')
	
	actual = db.get_keys()

	expected = ['id', 'startTime', 'endTime', 'status', 'a', 'b', 'c']

	assert actual == expected, "Cannot get update a row in the delimited file."

# Test if a cell in the delimited file can be updated
def test_DF_update_row():
	create_test_file()

	db = dapt.delimited_file.Delimited_file('test.csv', ',')
	
	db.update_row(1, {'id':'t2', 'startTime':'2019-09-06 17:37', 'endTime':'2019-09-06 17:55', 'status':'finished', 'a':'10', 'b':'10', 'c':'20'})
	
	expected = []
	expected.append(OrderedDict({'id':'t1', 'startTime':'2019-09-06 17:23', 'endTime':'2019-09-06 17:36', 'status':'finished', 'a':'2', 'b':'4', 'c':'6'}))
	expected.append({'id':'t2', 'startTime':'2019-09-06 17:37', 'endTime':'2019-09-06 17:55', 'status':'finished', 'a':'10', 'b':'10', 'c':'20'})
	expected.append({'id':'t3', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'-10', 'c':''})

	assert db.get_table() == expected, "Cannot update a row in the delimited file."

# Test if the row of a delimited file can be updated
def test_DF_update_cell():
	create_test_file()

	db = dapt.delimited_file.Delimited_file('test.csv', ',')
	
	db.update_cell(1, 'status', 'adding')
	
	expected = []
	expected.append(OrderedDict({'id':'t1', 'startTime':'2019-09-06 17:23', 'endTime':'2019-09-06 17:36', 'status':'finished', 'a':'2', 'b':'4', 'c':'6'}))
	expected.append({'id':'t2', 'startTime':'', 'endTime':'', 'status':'adding', 'a':'10', 'b':'10', 'c':''})
	expected.append({'id':'t3', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'-10', 'c':''})

	assert db.get_table() == expected, "Cannot update a cell in the delimited file."

# Test if the row of a delimited file can be updated
def test_DF_get_row_index():
	create_test_file()

	db = dapt.delimited_file.Delimited_file('test.csv', ',')

	assert db.get_row_index('status', 'finished') == 0, "Cannot get the row index."


