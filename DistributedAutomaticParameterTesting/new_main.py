'''
    Ben Duggan
    2/1/19
    Main script to run distributed parameter testing
'''

import dap
import os, sys


def main():
    conf = dap.config.Config('config.txt')
    sheet = dap.sheet.Sheet(conf.config['spreedsheetID'])
    ap = dap.param.Param(conf, sheet)
    boxy = dap.box.Box(conf)
    boxy.connect()

    #boxy.startServer() # Start the Flask server to get box access token

    print("Starting main script")

    while True:
        # Check the sheet for any trials that didn't run successfully
        ap.checkForDBErrors()

        parameters = ap.requestParameters() #Get the next parameter
        if parameters == None:
            print("No more parameters to run!")
            break

        print("Request parameters: ")
        print(parameters)

        file_name = 'dap_test_' + str(parameters['id']) + '.txt'

        f = open(file_name, 'w+')
        f.write(str(parameters))
        f.close()

        boxy.uploadFile(conf.config['boxFolderID'], '/', file_name)

        ap.parameterSuccessful(parameters["id"])

        exit()

        try:
            pass
            if 'clean' in parameters['tasks']:
                # Reset from the previous run
                print("Cleaning up folder")
                dap.tools.dataCleanup(conf.config)
                ap.updateStatus(parameters['id'], 'clean')

            if 'xml' in parameters['tasks']:
                # Create the parameters
                print("Creating parameters xml and autoParamSettings.txt")
                dap.tools.createXML(parameters)
                dap.tools.createSettingsFile(parameters)
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
                os.system('magick mogrify -format png *.svg')
                os.system('ffmpeg -framerate 24 -i snapshot%08d.png -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" output.mp4')
                os.chdir('../')

            if 'zip' in parameters['tasks']:
                # Zip Run output
                print("Zipping SVG and outputs")
                fileName = dap.tools.createZip(parameters)
                ap.updateStatus(parameters['id'], 'zip')

            if 'upload' in parameters['tasks']:
                # Upload zip to box
                if platform.system() == 'Windows':
                    uploadFile(ap.config['boxFolderID'], '\\', fileName)
                else:
                    uploadFile(ap.config['boxFolderID'], '/', fileName)

                ap.updateStatus(parameters['id'], 'upload')

            # Update sheets to mark the test is finished
            ap.parameterSuccessful(parameters["id"]) #Test completed successfully so we need to mark it as such

            # End tests
        except ValueError:
            # Test failed
            print(ValueError)
            print("Test failed")
            ap.parameterFailed(parameters["id"])

if __name__ == '__main__':
    #os.chdir("../")
    #print("Current working directory: ", os.getcwd())

    # Look at command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == 'reset':
            # Reset config.txt
            print("Reseting config file...")
            ap = autoParam()
            ap.config['userName'] = 'default'
            ap.config['numOfRuns'] = '1'
            ap.config['lastTest'] = 'None'
            ap.changeConfig()
            exit()

    main()
