"""
    Test if config.py is working correctly
"""

import dapt
import os

# Test creation of new config file
def test_config_create():
	dapt.config.Config.create('config.json')

	conf = dapt.config.Config('config.json')
	expected = {"last-test":None, "user-name":None, "spreedsheet-id":None, "client-id":None, "client-secret":None, "box-folder-id":None, "reset-time":None, "num-of-runs":None, "computer-strength":None, "access-token":None, "refresh-token":None}

	os.remove('config.json')

	assert conf.config == expected, "The Config.create() method that creates a new and empty config file returned wrong."

# Test changing a value in Config
def test_config_file_change():
	dapt.config.Config.create('config.json')

	conf = dapt.config.Config('config.json')

	expected = conf.config
	expected["user-name"] = 'Clifford'

	conf.config["user-name"] = 'Clifford'
	conf.update_config()
	conf.read_config()

	os.remove('config.json')

	assert conf.config == expected, "Config did not change the value correctly."

# Test adding a key,value pair in Config
def test_config_file_add():
	dapt.config.Config.create('config.json')

	conf = dapt.config.Config('config.json')

	expected = conf.config
	expected["abc"] = 123

	conf.config["abc"] = 123
	conf.update_config()
	conf.read_config()

	os.remove('config.json')

	assert conf.config == expected, "Config did not add the key,value pair correctly."

# Test making config file safe for uploading publicly
def test_config_safe():
	dapt.config.Config.create('config.json')
	conf = dapt.config.Config('config.json')

	expected = conf.config
	expected["access-token"] = ''
	expected["refresh-token"] = ''

	dapt.config.Config.safe('config.json')

	os.remove('config.json')

	assert conf.config == expected, "Conf did not make the config file safe."
