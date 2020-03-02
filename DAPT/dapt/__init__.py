

name = "dapt"
__version__ = "0.9.0.3"

import sys, argparse

from .config import Config
from .database import Database
from .delimited_file import Delimited_file
from .sheets import Sheet
from .box import Box
from .param import Param
from .tools import *

# Not used
def parse():
    parser = argparse.ArgumentParser(description='Distributed Automated Parameter Testing (DAPT)\nA library to assist with running parameter sets across multiple systems.', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--f', metavar='config.json', default='config.json', type=str, action='store', help="The path to the config file.")
    parser.add_argument('--r', action='store_true', help="Reset the config file.  \'last-test\':None")
    parser.add_argument('--c', action='store_true', help="Create a blank config file.")
    parser.add_argument('--s', action='store_true', help="Remove keys from the config file so it can be made public.")

    args = parser.parse_args()
    if args.r:
        # Remove last-test from config file
        conf = Config(args.f)
        if conf.config['last-test']:
            conf.config['last-test'] = None
        #if conf.config['performed-by']:
        #    conf.config['performed-by'] = None
        conf.update_config()
        exit()
    if args.c:
        # Reset config file
        #Config.create(args.f)
        exit()
    if args.s:
        # Safe config file
        Config.safe(args.f)
        exit()