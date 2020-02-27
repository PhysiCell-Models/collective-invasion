"""
Config
====== 

Class that allows for reading and modification of a configuration (config) file.  A config file is not required but using one will make using DAPT much easier to use and greatly increase increase it's functionality.  A configuration file is simply a JSON file.  There are some reserved keys but you can add your own and refer to them throughout your program.


Fields
^^^^^^

+---------------------------------+-----------------------------------------------------------------------------------------+
| Fields                          | Description                                                                             |
+=================================+=========================================================================================+
| ``last-test`` (str)             | The last test id that was run.                                                          |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``user-name`` (str)             | The box username of the user.                                                           |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``sheet-spreedsheet-id`` (str)  | The Google spreedsheet ID being used.                                                   |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``sheets-creds-path`` (str)     | The Google Sheets credentials file path.                                                |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``sheet-worksheet-id`` (str)    | The Google Sheets worksheet id.  Sheets are indexed at 0.                               |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``sheet-worksheet-title`` (str) | The Google Sheets worksheet title.                                                      |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``client-id`` (str)             | Box API client ID.                                                                      |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``client-secret`` (str)         | Box API client secret.                                                                  |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``box-folder-id`` (str)         | The box folder id to use                                                                |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``reset-time`` (str)            | The time that the box access-token needs to be refreshed.                               |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``num-of-runs`` (int)           | The number of paramater sets to run.                                                    |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``computer-strength`` (int)     | Any comments such as error messages relating to the parameter set.                      |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``access-token`` (str)          | The box access token for the particular session.                                        |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``refresh-token`` (str)         | The box refresh token for the particular session.                                       |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``remove-zip`` (bool)           | Have tools.data_cleanup() remove zip files if true.                                     |
+---------------------------------+-----------------------------------------------------------------------------------------+
| ``remove-movie`` (bool)         | Have tools.data_cleanup() remove mp4 files if true.                                     |
+---------------------------------+-----------------------------------------------------------------------------------------+


"""

import json

class Config:
    """
        Class which loads and allows for editing of a config file

        Args:
            path (string): path to config file
    """
    def __init__(self, path='config.json'):
        self.path = path
        self.config = self.read_config()

    def get_value(self, key, recursive=False):
        """
            Get the first value of the given key or return ``None`` if one doesn't exist.

            Args:
                key (str or list): the key (given as a string) or List containing the path to the value
                recursive (bool): recursively look through the config for the given key.  False by default.  If recursive is set to True then key must be a string.

            Returns:
                The value associated to the given key or None if the key is not in the dictionary.
        """

        # If we we want to do it recursively
        if recursive:
            return self._find_value(self.config, key)
        
        value = None

        # Convert key to list if it is not already one
        if not isinstance(key, list):
            key = [key]

        if key[0] in self.config:
            value = self.config[key[0]]

            for i in range(1, len(key)):
                if key[i] in value:
                    value = value[key[i]]
                else:
                    return None
        
        return value

    def _find_value(self, dic, key):
        """
            Helper function to find value recursively

            Args:
                dic (dict): the dictionary
                key (str): the key to search for

            Returns:
                The object associated to the given key or None if not found.
        """

        if key in dic:
            return dic[key]
        
        for k, v in dic.items():
            if isinstance(v, dict):
                value = self._find_value(v, key)
                if value:
                    return value
        return None

    
    def read_config(self):
        """
            Reads the file with path set to self.path

            Returns:
                Dictionary of config file
        """

        with open(self.path) as f:
            return json.load(f)

    def update_config(self, key=None, value=None):
        """
            Given a key and associated value, updated the config file.  Alternatively, you can give no arguments and the config dict will be saved.  You can also do both.

            Args:
                key (anything): the key to use.  If none is given then nothing will be updated in the dictionary.
                value (anything): the value associated ot the key.

            Returns:
                Dictionary of config file
        """

        if key:
            self.config[key] = value

        with open(self.path, 'w') as f:
            json.dump(self.config, f)

        return self.config

    @staticmethod
    def create(path='config.json'):
        """
            Creates a config file with the reserved keys inserted.

            Args:
                path (string): path where config file will be written
        """

        default = {"last-test":None, "user-name":None, "spreedsheet-id":None, "sheets-creds-path":None, "sheet-worksheet-id":None, "sheet-worksheet-title":None, "client-id":None, "client-secret":None, "box-folder-id":None, "reset-time":None, "num-of-runs":None, "computer-strength":None, "access-token":None, "refresh-token":None}
        
        with open(path, 'w') as f:
            json.dump(default, f)

    @staticmethod
    def safe(path="config.json"):
        """
            Safe config file by removing accessToken and refreshToken.

            Args:
                path (string): path where config file will be writen
        """
        conf = Config(path)
        data = conf.config
        if data["access-token"]:
            data["access-token"] = ""
        if data["refresh-token"]:
            data["refresh-token"] = ""
        conf.config = data
        conf.update_config()

