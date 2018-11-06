from __future__ import print_function
import os
import shutil
#import matlab.engine

def reprintline(line):
    print('\r', line, end='')

def runTests(tests):
    # Copy matlab scripts into output
    for file_name in os.listdir("matlab"):
        shutil.copy(os.getcwd()+"\\matlab\\"+file_name, "output")

if __name__ == '__main__':
    os.chdir("../")
    print("main")

    #runTests([])