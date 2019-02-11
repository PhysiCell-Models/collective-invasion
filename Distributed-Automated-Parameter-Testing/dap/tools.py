
import xml.etree.ElementTree as ET
import sys, os, platform, zipfile, datetime, time



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