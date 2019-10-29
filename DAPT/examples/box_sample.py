"""
    Example of how to use DAPT to upload results to Box.

    Before you can run this example, you need to create a new Box folder and record the paste the ID of that folder below.  If you are still confused you can follow the documentation guide.
"""

import dapt
import os, sys

# The box folder ID found in the URL of the folder
box_folder_id = '88499027100'

# These should be kept a secret!
client_id = 'ut798hvpbkj19pezrkgaptsk9dfxp6rf'
client_secret = 'PzwYQXLdtYYFIJNNTtR5TnuJItxOcveO'

# Create a Box object.  
# For this example we will directly pass the client ID and secret, however, these can be included in a config file and passed to Box via a Config object.
box = dapt.box.Box(client_id=client_id, client_secret=client_secret)

# Let's create a test file to upload to box
f = open('test.txt', 'w')
f.write('Box works correctly!!!')
f.close()

# In order to use box, a user must log in to get an access and refresh token.  
# You will need to go to 'http://127.0.0.1:5000' and login with your Box username and password to continue.
box.connect()

# Now we can upload the file to box
box.uploadFile(box_folder_id, '/test.txt', 'Box_test.txt')

# Lets clean up and remove the test file.
os.remove('test.txt')

