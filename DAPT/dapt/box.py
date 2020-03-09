"""
Box
=== 

Class that allows for access to the box API and methods to directly upload files
"""

import os, time
from boxsdk import *
from flask import *

class Box:
    """
    Class which allows for connection to box API.  You must either provide a Config object or client_id and client_secret.

    Keyword Args:
        config (Config): A Config object which contains the client_id and client_secret. 
        client_id (str): The Box client ID.
        client_secret (str): The Box client secret.
    """

    def __init__(self, *args, **kwargs):
        self.access_token = None
        self.refresh_token = None
        self.client = None
        self.refreshTime = None
        self.csrf_token = None

        self.app = Flask(__name__)

        if len(kwargs) == 0:
            raise ValueError("You must provide a Config object or client_id and client_secret.")
        elif len(kwargs) == 1:
            self.config=kwargs['config']

            self.oauth = OAuth2(
                client_id=self.config.config['client-id'],
                client_secret=self.config.config['client-secret']
            )
        elif len(kwargs) == 2:
            self.config = None

            self.oauth = OAuth2(
                client_id=kwargs['client-id'],
                client_secret=kwargs['client-secret']
            )
            
        else:
            raise ValueError("You entered too many arguments.")


    def connect(self, access_token = None, refresh_token = None):
        """
        Tries to connect to box using arguments provided in Config and starts server for authorization if not.

        Args:
            access_token (str): Optional argument that allows DAPT to connect to box without going through web authentification (assuming refresh_token is given and not expired).
            refresh_token (str): Optional argument that allows DAPT to connect to box without going through web authentification (assuming access_token is given and not expired).
        
        Returns:
            Box client if successful
        """

        # First check to see if the user gave us the access and refresh token
        if access_token and refresh_token:
            self.oauth._refresh_token = refresh_token
            self.access_token, self.refresh_token = self.oauth._refresh(access_token)
            self.client = Client(self.oauth)
            self.refreshTime = time.time() + 60*60

            return

        # If not, then we check to see if the access and refresh token are in the config file.
        if self.config and self.config.config['access-token'] and self.config.config['refresh-token']:
            try:
                print('Trying to get new access and refresh token from ' + self.config.path)
                self.oauth._refresh_token = self.config.config['refresh-token']
                self.access_token, self.refresh_token = self.oauth._refresh(self.config.config['access-token'])
                self.client = Client(self.oauth)
                self.refreshTime = time.time() + 60*60

                if self.config:
                    self.config.config['performed-by'] = self.client.user(user_id='me').get()['login'] #Save username to config
                    self.config.config['access-token'] = self.access_token #Save access token to config
                    self.config.config['refresh-token'] = self.refresh_token #Save refresn token to config
                    self.config.update_config()

                print('Got new access and refresh token from existing')
                return

            except Exception as e:
                print(e)

        print('Starting server.  Go to the URL below to activate box functionality.  If you are on a server you will need to run this code on your computer, get the access and refresh token and then add them to the config file.')
        return self.__startServer()

    def __startServer(self):
        """
        Method that starts flask to start authorization process

        Returns:
            Box client which can be used to access authorized user data
        """

        print("Starting server.  Go to 127.0.0.1:5000 to authenticate box.  It can only be ended by completing authentification or going to 127.0.0.1:5000/end")
        self.app.add_url_rule('/', 'index', self.__index)
        self.app.add_url_rule('/return', 'return', self.__capture)
        self.app.add_url_rule('/end', 'end', self.__end)
        self.app.run()
        print("Server stoped")
        return self.client

    def __index(self):
        """
        Flask page: index of the web server and serves as the start point for authentication

        Returns:
            String containing HTML to be displayed
        """
        self.auth_url, self.csrf_token = self.oauth.get_authorization_url("http://127.0.0.1:5000/return")

        return '<h1>Welcome to box auth</h1> This web server is used to interface with the box API.  Click the link below to securely login on box.' + '<a href="'+self.auth_url+'">Click here to authenticate your box account </a>'
    
    def __capture(self):
        """
        Flask page: box redirect url which contains the code and state used to get access and refresh token

        Returns:
            String containing HTML to be displayed with box login credentials
        """

        # Capture auth code and csrf token via state
        code = request.args.get('code')
        state = request.args.get('state')

        # If csrf token matches, fetch tokens
        print(self.csrf_token)
        assert state == self.csrf_token
        self.access_token, self.refresh_token = self.oauth.authenticate(code)

        self.client = Client(self.oauth)

        self.refreshTime = time.time() + 60*60

        func = request.environ.get('werkzeug.server.shutdown')
        if func is None:
            raise RuntimeError('Not running with the Werkzeug Server')
        func()

        if self.config:
            self.config.config['performed-by'] = self.client.user(user_id='me').get()['login'] #Save username to config
            self.config.config['access-token'] = self.access_token #Save access token to config
            self.config.config['refresh-token'] = self.refresh_token #Save refresn token to config
            self.config.update_config()


        return 'You are now logged in as: ' + self.client.user(user_id='me').get()['login'] + '<br><strong>The server has been shutdown and the normal script is resuming.</strong><br>access token: '+self.access_token+'<br>refresh token: '+self.refresh_token+'<br><a href="http://127.0.0.1:5000">Click to go to index (assuming server restarted)</a>'

    def __end(self):
        """
        Flask page: shuts down flask server

        Returns:
            String containing HTML to be displayed
        """

        func = request.environ.get('werkzeug.server.shutdown')
        if func is None:
            raise RuntimeError('Not running with the Werkzeug Server')
        func()
        return 'Server shutting down...<br> Getting back to the main python script.'

    def updateTokens(self, access_token):
        """
        Refresh the access and refresh token given a valid access token

        Args:
            access_token (string): box access token to be refreshed

        Returns:
            Box client
        """

        self.access_token, self.refresh_token = self.oauth.refresh(access_token)
        self.client = Client(self.oauth)
        self.refreshTime = time.time() + 60*60

        if self.config:
            self.config.config['performed-by'] = self.client.user(user_id='me').get()['login'] #Save username to config
            self.config.config['access-token'] = self.access_token #Save access token to config
            self.config.config['refresh-token'] = self.refresh_token #Save refresn token to config
            self.config.update_config()

        return self.client

    def uploadFile(self, folderID, path, name):
        """
        Upload a file to box using the current client

        Args:
            folderID (string): the box folder ID that the file should be added to (found in URL)
            path (string): the path to the file on your local machine
            file (string): the name of the file you wish to upload
        """

        # If the access token is expired (or about to be), then update it
        if self.refreshTime < time.time() + 5:
            self.updateTokens(self.access_token)
        print(os.getcwd()+path, name)
        return self.client.folder(folderID).upload(os.getcwd()+path, name)

