.. _install:

Install
=======

To install DAPT you can either use a `release <https://github.com/BenSDuggan/DAPT/releases>`_ downloaded from GitHub (recomended), or download the most recent "stable" version of the code by cloning or downloading from the master brantch.  If you install using a release they can be found at `https://github.com/BenSDuggan/DAPT/releases <https://github.com/BenSDuggan/DAPT/releases>`_.  You can also download from the GitHub page for run ``git clone https://github.com/BenSDuggan/DAPT``.

Once you've downloaded the code, you can install the python libraries by running ``pip3 install -r requirements.txt`` in the terminal from inside the root folder. You can then open a python session and test to see if everything installed correctly:

.. code-block:: python
	
	import dapt
	dapt.__version__

You should then see some version number corresponding with the version you downloaded.


You can now use the library!  However, the functionality can be greatly increased by connecting some other services such as Google Sheets or Box.  See the below guides on how to include this functionality.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   sheets-install
   box-install