"""
    Config
    ====== 

    Class that allows for reading and modification of a configuration (config) file.  A config file is not required but using one will make using DAPT much easier to use and greatly increase increase it's functionality.  A configuration file is simply a JSON file.  There are some reserved keys but you can add your own and refer to them throughout your program.


    Fields
    ^^^^^^

    +-----------------------------+-----------------------------------------------------------------------------------------+
    | Fields                      | Description                                                                             |
    +=============================+=========================================================================================+
    | ``last-test`` (str)         | The last test id that was run.                                                          |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``user-name`` (str)         | The box username of the user.                                                           |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``spreedsheet-id`` (str)    | The Google spreedsheet ID being used.                                                   |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``client_id`` (str)         | Box API client ID.                                                                      |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``client_secret`` (str)     | Box API client secret.                                                                  |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``box-folder-id`` (str)     | The box folder id to use                                                                |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``reset-time`` (str)        | The time that the box access-token needs to be refreshed.                               |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``num-of-runs`` (int)       | The number of paramater sets to run.                                                    |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``computer-strength`` (int) | Any comments such as error messages relating to the parameter set.                      |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``access-token`` (str)      | The box access token for the particular session.                                        |
    +-----------------------------+-----------------------------------------------------------------------------------------+
    | ``refresh-token`` (str)     | The box refresh token for the particular session.                                       |
    +-----------------------------+-----------------------------------------------------------------------------------------+


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
    
    def read_config(self):
        """
            Reads the file with path set to self.path

            Returns:
                Dictionary of config file
        """

        with open(self.path) as f:
            return json.load(f)

    def update_config(self):
        """
            Change value for a given key in the config file

            Returns:
                Dictionary of config file
        """

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

        default = {"last-test":None, "user-name":None, "spreedsheet-id":None, "client-id":None, "client-secret":None, "box-folder-id":None, "reset-time":None, "num-of-runs":None, "computer-strength":None, "access-token":None, "refresh-token":None}
        
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

