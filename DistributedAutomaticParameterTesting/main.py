# Ben Duggan
# 10/15/18
# Main script to run distributed parameter testing

from __future__ import print_function
import xml.etree.ElementTree as ET
from autoParam import *
import os

def reprintline(line):
    print('\r', line, end='')

def createXML(parameters, offLimits=[]):
    parameters = dict(parameters)
    tree = ET.parse("config/PhysiCell_settings_default.xml")
    root = tree.getroot()

    for key in parameters:
        if key in offLimits:
            del parameters[key]

    for key in parameters:
        node = root.find(key)

        if node != None:
            node.text = str(parameters[key])

    tree.write("config/PhysiCell_settings.xml")

def createSettingsFile(parameters):
    data = ""

    for key in parameters:
        data += str(key) + ":" + str(parameters[key]) + "\n"

    with open('autoParamSettings.txt', 'w') as file:
        file.writelines(data)

def main():
    ap = autoParam()
    count = 0

    while count < int(ap.config["numOfRuns"]) or int(ap.config["numOfRuns"]) == -1:
        # Check the sheet for any trials that didn't run successfully
        ap.checkForDBErrors()

        parameters = ap.requestParameters() #Get the next parameter
        if parameters == None:
            print("No more parameters to run!")
            break

        print("Request parameters: ")
        print(parameters)

        try:
            # Start tests
            createXML(parameters)
            createSettingsFile(parameters)


            # End tests
        except:
            # Test failed
            print("Test failed")
            ap.parameterFailed(parameters["id"])

        # Update sheets to mark the test is finished
        ap.parameterSuccessful(parameters["id"]) #Test completed successfully so we need to mark it as such

        count += 1

if __name__ == '__main__':
    os.chdir("../")
    print("Current working directory: ", os.getcwd())

    main()

    '''
    ap = autoParam()
    count = 0

    parameters = ap.requestParameters()
    createXML(parameters)
    createSettingsFile(parameters)
    '''
