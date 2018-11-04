'''
    Ben Duggan
    10/30/18
    Main script to run distributed parameter testing
'''

import xml.etree.ElementTree as ET
from autoParam import *
import os
import platform
import zipfile
import datetime
import time
#import matlab.engine
from box import *

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

def dataCleanup(config):
    # Emulating make data-cleanup: remove .mat, .xml, .svg, .txt, .pov
    for file in os.listdir("."):
        if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".pov") or (config['removeZip']=='True' and file.endswith('.zip')):
            os.remove(file)

    for file in os.listdir("SVG/"):
        #if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".pov"):
        os.remove("SVG/" + file)

    for file in os.listdir("output/"):
        #if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".pov"):
         os.remove("output/" + file)

def createSettingsFile(parameters):
    data = ""

    for key in parameters:
        data += str(key) + ":" + str(parameters[key]) + "\n"

    with open('autoParamSettings.txt', 'w') as file:
        file.writelines(data)

def createZip(parameters):
    fileName = str(parameters['id']) + '_test_' + datetime.datetime.utcnow().strftime('%Y-%m-%d_%H-%M-%S') + '.zip'
    # Create the zip
    zip = zipfile.ZipFile(fileName, 'w')

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

    return fileName


def main():
    ap = autoParam()
    count = 0

    useBox = None
    while(useBox != "Y" and useBox != "y" and useBox != "n" and useBox != "N"):
        useBox = input("Do you want to use box (y/n): ")
    if useBox == "Y" or useBox == "y":
        useBox = True
        username = startServer()
        ap.config['userName'] = username
        ap.changeConfig()
    print("Starting main script")

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
            if 'clean' in parameters['tasks']:
                # Reset from the previous run
                print("Cleaning up folder")
                dataCleanup(ap.config)
                ap.updateStatus(parameters['id'], 'clean')

            if 'xml' in parameters['tasks']:
                # Create the parameters
                print("Creating parameters xml and autoParamSettings.txt")
                createXML(parameters)
                createSettingsFile(parameters)
                ap.updateStatus(parameters['id'], 'xml')

            if 'sim' in parameters['tasks']:
                # Run PhysiCell
                print("Running test")
                if platform.system() == 'Windows':
                    os.system("AMIGOS-invasion.exe")
                else:
                    os.system("./AMIGOS-invasion")
                ap.updateStatus(parameters['id'], 'sim')

            if 'matlab' in parameters['tasks']:
                # Run matlab scripts
                print("Run matlab scripts")
                print("^^ not implimented ^^")

            if 'imgProc' in parameters['tasks']:
                # Run image processing
                print("Run image processing")
                os.chdir('output/')
                os.system('mogrify -format png *.svg')
                os.system('ffmpeg -framerate 24 -i snapshot%08d.png -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" output.mp4')
                os.chdir('../')

            if 'zip' in parameters['tasks']:
                # Zip Run output
                print("Zipping SVG and outputs")
                fileName = createZip(parameters)
                ap.updateStatus(parameters['id'], 'zip')

            if 'upload' in parameters['tasks']:
                # Upload zip to box
                if useBox:
                    print("Uploading zip to box")
                    if platform.system() == 'Windows':
                        uploadFile(ap.config['boxFolderID'], '\\', fileName)
                    else:
                        uploadFile(ap.config['boxFolderID'], '/', fileName)

                    ap.updateStatus(parameters['id'], 'upload')
                else:
                    print("Cannot upload to box")

            # End tests
        except ValueError:
            # Test failed
            print(ValueError)
            print("Test failed")
            ap.parameterFailed(parameters["id"])

        # Update sheets to mark the test is finished
        ap.parameterSuccessful(parameters["id"]) #Test completed successfully so we need to mark it as such

        count += 1

if __name__ == '__main__':
    os.chdir("../")
    print("Current working directory: ", os.getcwd())

    main()
