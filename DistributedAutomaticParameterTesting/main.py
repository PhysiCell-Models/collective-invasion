# Ben Duggan
# 10/21/18
# Main script to run distributed parameter testing

from __future__ import print_function
import xml.etree.ElementTree as ET
from autoParam import *
import os
import zipfile
import datetime

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

def dataCleanup():
    '''
    Emulating make data-cleanup
        rm -f *.mat
        rm -f *.xml
        rm -f *.svg
        rm -f *.txt
        rm -f *.pov
        rm -f ./Output/*
        rm -f ./SVG/*
    '''

    # remove .mat, .xml, .svg, .txt, .pov
    for file in os.listdir("."):
        if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".pov"):
            os.remove(file)

    for file in os.listdir("SVG/"):
        if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".pov"):
            os.remove("SVG/" + file)

    for file in os.listdir("output/"):
        if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".pov"):
            os.remove("output/" + file)


def createSettingsFile(parameters):
    data = ""

    for key in parameters:
        data += str(key) + ":" + str(parameters[key]) + "\n"

    with open('autoParamSettings.txt', 'w') as file:
        file.writelines(data)

def createZip(parameters):
    # Create the zip
    zip = zipfile.ZipFile(str(parameters['id']) + '_test_' + datetime.datetime.utcnow().strftime('%Y-%m-%d_%H-%M-%S') + '.zip', 'w')

    # Add individual files
    zip.write('autoParamSettings.txt', compress_type=zipfile.ZIP_DEFLATED)

    # Add files in SVG
    for folder, subfolders, files in os.walk('SVG/'):
        for file in files:
            zip.write(os.path.join(folder, file), 'SVG/'+file, compress_type = zipfile.ZIP_DEFLATED)

    # Add files in output
    for folder, subfolders, files in os.walk('output/'):
        for file in files:
            zip.write(os.path.join(folder, file), 'output/'+file, compress_type = zipfile.ZIP_DEFLATED)

    zip.close()

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
            # Reset from the previous run
            print("Cleaning up folder")
            dataCleanup()

            # Create the parameters
            print("Creating parameters xml and autoParamSettings.txt")
            createXML(parameters)
            createSettingsFile(parameters)

            # Run PhysiCell
            print("Running test")
            os.system("AMIGOS-invasion.exe")

            # Zip Run output
            print("Zipping SVG and outputs")
            createZip(parameters)

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
