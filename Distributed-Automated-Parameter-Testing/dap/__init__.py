

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
    exit()

from . import sheet
from . import param
from . import box
from . import tools
