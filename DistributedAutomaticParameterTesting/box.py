'''
    Ben Duggan
    10/30/18
    Script for getting box access token and uploading files
'''

import os
from boxsdk import *
from flask import *
import time

''' START - Flask web server used for auth'''
app = Flask(__name__)
global access_token, refresh_token, client, refreshTime

# The method called to start the server
def startServer():
    print("Starting server.  Go to 127.0.0.1:5000 to authenticate box.  It can only be ended by completing authentification or going to 127.0.0.1:5000/end")
    app.run()
    print("Server stoped")
    return client.user(user_id='me').get()['login']

# This is the index of the web server and serves as the start point for authentication
@app.route('/')
def index():
    global oauth
    oauth = OAuth2(
        client_id=config['client_id'],
        client_secret=config['client_secret']
    )

    global csrf_token
    auth_url, csrf_token = oauth.get_authorization_url("http://127.0.0.1:5000/return")

    return '<h1>Welcome to box auth</h1> This web server is used to interface with the box API.  Click the link below to securely login on box.' + '<a href="'+auth_url+'">Click here to authenticate your box account </a>'

# This page is called using a GET request by the Box API and gives the access token and refresh token needed to upload files to box
@app.route('/return')
def capture():
    # Capture auth code and csrf token via state
    code = request.args.get('code')
    state = request.args.get('state')

    # If csrf token matches, fetch tokens
    assert state == csrf_token
    global access_token
    global refresh_token
    access_token, refresh_token = oauth.authenticate(code)

    global client
    client = Client(oauth)

    global refreshTime
    refreshTime = time.time() + 60*60

    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()

    return 'You are now logged in as: ' + client.user(user_id='me').get()['login'] + '<br><strong>The server has been shutdown and the normal script is resuming.</strong><a href="http://127.0.0.1:5000">Click to go to index (assuming server restarted)</a>'

# When the user goes to /end the server shuts down
@app.route('/end')
def end():
    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()
    return 'Server shutting down...<br> Getting back to the main python script.'
''' END - Flask web server used for auth'''

''' START - Box upload functions '''
f = open(''+'config.txt', 'r').readlines()
global config
config = {}

for i in range(0, len(f)):
    f[i] = f[i].replace("\n", "")
    config[f[i].split(":")[0]] = f[i].split(":")[1]

print(config)


def uploadFile(folderID, path, file):
    if refreshTime < time.time() + 5:
        updateTokens(access_token)
    print(client.folder(folderID).upload(os.getcwd()+path+file, file))

def updateTokens(access):
    global access_token, refresh_token, client, refreshTime
    access_token, refresh_token = oauth.refresh(access)
    client = Client(oauth)
    refreshTime = time.time() + 60*60

def getAccess():
    return access_token
def getRefresh():
    return refresh_token
def getClient():
    return client

''' END - Box upload functions '''

if __name__ == '__main__':
    os.chdir("../")
    app.run()

    items = client.folder(folder_id='0').get_items(limit=100, offset=0)

    for i in items:
        print(i)


    #uploadFile('\\DistributedAutomaticParameterTesting\\', "testPayload.zip")