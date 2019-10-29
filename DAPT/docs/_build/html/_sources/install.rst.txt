.. _install:

Install
=======

DAPT
----

To install DAPT you can either use a `release <https://github.com/BenSDuggan/DAPT/releases>`_ downloaded from GitHub (recomended) or download the most recent version of the code.  If you install using a release they can be found at `https://github.com/BenSDuggan/DAPT/releases <https://github.com/BenSDuggan/DAPT/releases>`_.  Once downloaded unzip the data.  You can also download from the GitHub page for run ``git clone https://github.com/BenSDuggan/DAPT``.

Once you've downloaded the code, you can install the python libraries by running ``pip install -r requirements.txt`` in the terminal from inside the root folder. You can then open a python session and test to see if everything installed correctly:

.. code-block:: python
	
	import dapt
	dapt.__version__


Google Sheets
-------------

Google Sheets can be used as a database to store your paramater sets.  The advantage to using this "database" over a file containing the paramaters is that a team can work from and update the paramater list.  

API Credentials
^^^^^^^^^^^^^^^

To use Google Sheets you will need to use the Google Sheets API and generate the proper credentials.  Start by going to the `Google Developer Console <https://console.developers.google.com>`_ and login using a Google account.  Next create a new project and name it whatever you'd like.  Next click "ENABLE APIS AND SERVICES", search for "Google Sheets API" and click it.  Then click "Enable".  Then click "Credentials" from the menu on the left side of the page.

Click the dropdown at the top of the page that says "CREATE CREDENTIALS" and select "Service account".  Give the service a name and click "Create".  In the next section asking about permissions select "Project" then "Editor" and then select "Continue".  You don't need to add anyone in the optional 3rd step.  Select "CREATE KEY", ensure the key type of "JSON" is selected and click "CREATE".  You will give DAPT the path to this JSON file when using Google Sheets.

Preparing a Sheet
^^^^^^^^^^^^^^^^^


Box
---

API Credentials
^^^^^^^^^^^^^^^

Preparing a folder
^^^^^^^^^^^^^^^^^^

