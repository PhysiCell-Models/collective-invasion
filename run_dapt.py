"""
Ben Duggan modified by John Metzcar for ECM invasion project
10/22/2020
Main script to run the DAPT pipeline to test parameters.

* To run this code, the configuration file (named config.json) and Google Sheet credentials (credentials.json) should be located one direcotry up.
* You can upload to Box or Google Drive by setting the boolean True/False in the `box_upload` and `google_drive_upload`, respectively.
"""

import os, platform, datetime, zipfile
import dapt

# Pipeline variables
box_upload = False
google_drive_upload = True
config_path = '../configs/config.json'
google_sheets_creds_path = '../configs/google-sheets-creds.json'
google_drive_creds_path = '../configs/google-drive-creds.json'

# Initialize DAPT
conf = dapt.Config(config_path)
sheet = dapt.Sheet(config=conf, creds=google_sheets_creds_path)
ap = dapt.Param(sheet, conf)

# Where should we upload
if box_upload:
    boxy = dapt.storage.Box(config=conf)
    boxy.connect(access_token=conf["box"]["access-token"], refresh_token=conf["box"]["refresh-token"])
if google_drive_upload:
    drive = dapt.storage.Google_Drive(creds_path=google_drive_creds_path, config=conf)
    drive.connect()

def createZip(parameters):
    fileName = str(parameters['id']) + '_test_' + datetime.datetime.utcnow().strftime('%Y-%m-%d_%H-%M-%S') + '.zip'
    # Create the zip
    zip = zipfile.ZipFile(fileName, 'w')

    # Add individual files
    zip.write(parameters["id"]+'_dapt_param_settings.txt', compress_type=zipfile.ZIP_DEFLATED)

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

print("Starting main script")

parameters = ap.next_parameters()

while parameters:
    print("Request parameters: ")
    print(parameters)

    try:
        if 'clean' in parameters['tasks']:
            # Reset from the previous run
            print("Cleaning up folder")
            ap.update_status(parameters['id'], 'clean')

            dapt.data_cleanup(config=conf)

        if 'xml' in parameters['tasks']:
            # Create the parameters
            print("Creating parameters xml and autoParamSettings.txt")
            ap.update_status(parameters['id'], 'xml')

            dapt.create_XML(parameters, default_settings="config/PhysiCell_settings_default.xml", save_settings="config/PhysiCell_settings.xml")
            dapt.create_settings_file(parameters)

        if 'sim' in parameters['tasks']:
            # Run PhysiCell
            print("Running test")
            ap.update_status(parameters['id'], 'sim')
            
            if platform.system() == 'Windows':
                os.system("AMIGOS-invasion.exe")
            else:
                os.system("./AMIGOS-invasion")

        if 'matlab' in parameters['tasks']:
            # Run matlab scripts
            print("Run matlab scripts")
            ap.update_status(parameters['id'], 'matlab')
            print("^^ not implimented ^^")

        if 'imgProc' in parameters['tasks']:
            # Run image processing
            print("Run image processing")
            ap.update_status(parameters['id'], 'imgProc')

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
            ap.update_status(parameters['id'], 'zip')

            zipName = createZip(parameters)
            

        if 'upload' in parameters['tasks']:
            # Upload to box or google drive
            print("Upload zip")
            ap.update_status(parameters['id'], 'upload')

            if box_upload:
                # Upload zip to box
                print("Uploading zip to box")

                print(boxy.upload_file(conf['boxFolderID'], fileName))
                print(boxy.upload_file(conf['boxFolderZipID'], zipName))
            if google_drive_upload:
                # Upload zip to Google Drive
                print("Uploading zip to Google Drive")

                #print(drive.upload_file(conf["drive-folder-id"], fileName))
                print(drive.upload_file(conf["drive-zip-id"], zipName))

        # Update sheets to mark the test is finished
        ap.successful(parameters["id"]) #Test completed successfully so we need to mark it as such
        
    except ValueError:
        # Test failed
        print("Test failed")
        print(ValueError)
        ap.failed(parameters["id"], ValueError)
    
    parameters = ap.next_parameters()

print("Finished running parameters!")