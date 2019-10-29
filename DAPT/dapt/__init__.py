"""

Distributed Automated Paramater Testing (DAPT)
==============================================

This is __init__.py
"""

__version__ = "0.8.2"

import sys
from . import config

# Parse command line arguments
if len(sys.argv) > 1:
    # Config settings
    if sys.argv[1] == 'config':
        file_name = 'c.txt'

        # Check for flags
        try:
            if '-f' in sys.argv[1:]:
                file_name = sys.argv[sys.argv.index('-f')+1]
        except:
            print('Couldn\'t use specified file.  Using file name "config.txt".')

        # Reset config file
        if sys.argv[2] == 'create':
            config.Config.create(file_name)

        # Safe config file
        if sys.argv[2] == 'safe':
            config.Config.safe(file_name)
    #exit()

from . import box
from . import database
from . import delimited_file
from . import sheets
from . import param
from . import tools

def test():
    from . import tests
    import pytest

    pytest.main()