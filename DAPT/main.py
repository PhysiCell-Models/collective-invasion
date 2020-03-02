'''
    Ben Duggan modified by John Metzcar for ECM invasion project
    06/24/19
    Main script to run distributed parameter testing
'''

import xml.etree.ElementTree as ET
import sys,os,platform,zipfile,datetime, time
import dapt

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
    # Emulating make data-cleanup-light: remove .mat, .xml, .svg, .txt, .pov
    for file in os.listdir("."):
        if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".pov") or file.endswith(".png") or (config['removeZip']=='True' and file.endswith('.zip')) or (config['removeMovie']=='True' and file.endswith('.mp4')):
            os.remove(file)

    for file in os.listdir("output/"):
        if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".png"):
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

    # Add files in output
    for folder, subfolders, files in os.walk('output/'):
        for file in files:
            if ".png" not in file:
                zip.write(os.path.join(folder, file), 'output/'+file, compress_type = zipfile.ZIP_DEFLATED)

    # Add programming files
    zip.write('config/PhysiCell_settings.xml', compress_type = zipfile.ZIP_DEFLATED)
    zip.write('main-ecm.cpp', compress_type = zipfile.ZIP_DEFLATED)
    zip.write('Makefile', compress_type = zipfile.ZIP_DEFLATED)
    zip.write('custom_modules/AMIGOS-invasion.h', compress_type = zipfile.ZIP_DEFLATED)
    zip.write('custom_modules/AMIGOS-invasion.cpp', compress_type = zipfile.ZIP_DEFLATED)

    zip.close()

    return fileName


def main():
    conf = dapt.Config('DAPT/config.json')
    sheet = dapt.Sheet(config=conf, creds='DAPT/credentials.json')
    ap = dapt.Param(sheet, conf)
    boxy = dapt.Box(config = conf)
    boxy.connect()

    print("Starting main script")

    while True:
        # Check the sheet for any trials that didn't run successfully
        #ap.checkForDBErrors()

        parameters = ap.next_parameters() #Get the next parameter
        if parameters == None:
            print("No more parameters to run!")
            break

        print("Request parameters: ")
        print(parameters)

        try:
            if 'clean' in parameters['tasks']:
                # Reset from the previous run
                print("Cleaning up folder")
                dataCleanup(conf.config)
                ap.update_status(parameters['id'], 'clean')

            if 'xml' in parameters['tasks']:
                # Create the parameters
                print("Creating parameters xml and autoParamSettings.txt")
                createXML(parameters)
                createSettingsFile(parameters)
                ap.update_status(parameters['id'], 'xml')

            if 'parameters' in parameters['tasks']:
                createSettingsFile(parameters)
                ap.update_status(parameters['id'], 'parameters')

            if 'sim' in parameters['tasks']:
                # Run PhysiCell
                print("Running test")
                if platform.system() == 'Windows':
                    os.system("AMIGOS-invasion.exe")
                else:
                    os.system("./AMIGOS-invasion")
                ap.update_status(parameters['id'], 'sim')

            if 'matlab' in parameters['tasks']:
                # Run matlab scripts
                print("Run matlab scripts")
                print("^^ not implimented ^^")

            if 'imgProc' in parameters['tasks']:
                # Run image processing
                print("Run image processing")
                fileName = str(parameters['id']) + '_test_' + datetime.datetime.utcnow().strftime('%Y-%m-%d_%H-%M-%S') + str('.mp4')
                print(fileName)
                os.chdir('output/')
                os.system('mogrify -format png *.svg')
                movie_run_command_str = str('ffmpeg -framerate 24 -i snapshot%08d.png -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" ../')
                movie_run_command_str = movie_run_command_str + fileName
                os.system(movie_run_command_str)
                os.chdir('../')

            if 'zip' in parameters['tasks']:
                # Zip Run output
                print("Zipping SVG and outputs")
                zipName = createZip(parameters)
                ap.update_status(parameters['id'], 'zip')

            if 'upload' in parameters['tasks']:
                # Upload zip to box
                print("Uploading zip to box")
                print(zipName)
                if platform.system() == 'Windows':
                    print(boxy.uploadFile(conf.config['boxFolderID'], '\\', fileName))
                    print(boxy.uploadFile(conf.config['boxFolderZipID'], '\\', zipName))
                else:
                    print(boxy.uploadFile(conf.config['boxFolderID'], '/'+fileName, fileName))
                    print(boxy.uploadFile(conf.config['boxFolderZipID'], '/'+zipName, zipName))

                ap.update_status(parameters['id'], 'upload')

            # Update sheets to mark the test is finished
            ap.successful(parameters["id"]) #Test completed successfully so we need to mark it as such
            
            # End tests
        except ValueError:
            # Test failed
            print(ValueError)
            print("Test failed")
            ap.failed(parameters["id"], ValueError)
        
if __name__ == '__main__':
    os.chdir("../")
    print("Current working directory: ", os.getcwd())

    main()