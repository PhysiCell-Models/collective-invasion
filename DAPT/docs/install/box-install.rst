.. _box-install:

Box Installation
================

Box is a cloud storage service that many universities allow students, facutie and staff to use.  The advantage of box is that it allows a large amount of data to be uploaded to a common place where team memebers can observe data.  In order to allow DAPT to upload to box, you must create some API credentials.


API Credentials
---------------

1. Start by going to the `Box Development <https://developer.box.com/>`_ website and clicking on the blue "Console" button.  Then log in.
2. Click "Create New App".  Then click "Custom App" and "Next" on the next page.
3. On the "Authentication Method" page click "Standard OAuth 2.0 (User Authentication)" and name your project.  Then click "View Your App".
4. Scroll down to the "OAuth 2.0 Credentials" section and record the Client ID and Secret.  You will pass these to the DAPT Box class to allow the Box SDK to work.
5. Lastely, scroll down to the "OAuth 2.0 Credentials" section and change the url to ``http://127.0.0.1:5000/return``.  Then click "Save Changes".

