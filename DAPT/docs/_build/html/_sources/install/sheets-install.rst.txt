.. _sheets-install:

Google Sheets Installation
==========================

Google Sheets can be used as a database to store your paramater sets.  The advantage to using this "database" over a file containing the paramaters is that a team can work on the set more collabortively and update the paramater list on the fly.  The bigger advantage is that the paramater list can be run dynamically.  By this I mean that people, running the library simultaneously, can connect to Google Sheets and get the next paramater set in the list.

API Credentials
---------------

To use Google Sheets you will need to use the Google Sheets API and generate the proper credentials.  

1. Start by going to the `Google Developer Console <https://console.developers.google.com>`_ and login using a Google account.  
2. Create a new project and name it whatever you'd like.  You can do this by selecting the down arrow next in the top left corner of the page next to the "Google API" logo and clicking "New Project".
3. Click "ENABLE APIS AND SERVICES", search for "Google Sheets API" and click it.  Then click "Enable".  
4. Click the "Credentials" tab from the menu on the left side of the page.
5. Click the dropdown at the top of the page that says "CREATE CREDENTIALS" and select "Service account".  
6. Give the service a name and click "Create".  
7. In the next section asking about service account permissions, create a role with by selecting "Project" then "Editor".  Then select "Continue".
8. Select "CREATE KEY", ensure the key type of "JSON" is selected and click "CREATE".  You will give DAPT the path to this JSON file when using Google Sheets.  Then click "DONE".
9. You should now be on the Credentials page of the Google Sheets API.  Record the email address in the Service Account table.  It should look something like ``*.iam.gserviceaccount.com``.
