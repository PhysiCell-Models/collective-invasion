import os
import zipfile
from datetime import datetime

def createZip(parameters):
    # Create the zip
    #datetime.utcnow().isoformat()
    zip = zipfile.ZipFile(str(parameters['id']) + '_test_' + "d" +'.zip', 'w')

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

if __name__ == '__main__':
    os.chdir("../")
    #os.system("AMIGOS-invasion.exe")
    #createZip({"id":8})
    print(datetime.utcnow().isoformat())
