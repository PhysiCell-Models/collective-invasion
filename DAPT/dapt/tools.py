"""
Tools
=====

A collection of tools that make DAPT easy to use with PhysiCell.
"""

import xml.etree.ElementTree as ET
import sys, os, platform, zipfile, datetime, time


def create_XML(parameters, default_settings="PhysiCell_settings_default.xml", save_settings="PhysiCell_settings.xml", off_limits=[]):
    """
        Create a PhysiCell XML settings file given a dictionary of paramaters.  This function works by having a ``default_settings`` file which contains the generic XML structure.  Each key in ``parameters` then contains the paths to each XML tag in the ``default_settings`` file.  The value of that tag is then set to the value in the associated key.  If a key in ``parameters`` does not exist in the ``default_settings`` XML file then it is ignored.  If a key in ``parameters`` also exists in ``off_limits`` then it is ignored.

        Args:
            paramaters (dict): A dictionary of paramaters where the key is the path to the xml variable and the value is the desired value in the XML file.
            default_settings (str): the path to the default xml file
            save_settings (str): the path to the output xml file
            off_limits (list): a list of keys that should not be inserted into the XML file.
    """

    parameters = dict(parameters)
    tree = ET.parse(default_settings)
    root = tree.getroot()

    for key in parameters:
        if key in off_limits:
            next

        node = root.find(key)

        if node != None:
            node.text = str(parameters[key])

    tree.write(save_settings)

def data_cleanup(config=None):
    """
        Emulating make data-cleanup-light: remove .mat, .xml, .svg, .txt, .pov.  You can optionally remove zipped files by setting ``remove-zip`` equal to ``True`` or remove ``*.mp4`` by setting ``remove-movie`` to ``True`` in the config file.

        Args:
            config (Config): A config object, optionally given.
    """

    for file in os.listdir("."):
        if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".pov") or file.endswith(".png"):
            os.remove(file)
        elif config:
            if (config.get_value('remove-zip', recursive=True) and file.endswith('.zip')) or (config.get_value('remove-movie', recursive=True) and file.endswith('.mp4')):
                os.remove(file)

    for file in os.listdir("output/"):
        if file.endswith(".mat") or file.endswith(".xml") or file.endswith(".svg") or file.endswith(".txt") or file.endswith(".png"):
            os.remove("output/" + file)
        elif config:
            if (config.get_value('remove-zip', recursive=True) and file.endswith('.zip')) or (config.get_value('remove-movie', recursive=True) and file.endswith('.mp4')):
                os.remove("output/" + file)

def create_settings_file(parameters, pid=None):
    """
        Creates a file where each line contains a key from the parameter file and the associated key, separated by a semicolon.

        Args:
            parameters (dict): the paramaters to be saved in the file
            pid (str): the parameter id of the current parameter run.  If you don't give an id then the id in ``parameters`` will be used.
    """

    data = ""

    if not pid:
        pid = parameters["id"]

    for key in parameters:
        data += str(key) + ":" + str(parameters[key]) + "\n"

    with open(str(parameters["id"]) + "_dapt_param_settings.txt", 'w') as file:
        file.writelines(data)

def create_zip(pid):
    """
        Zip all of the important PhysiCell items.

        Args:
            pid (str): the id of the current parameter run

        Returns:
            The name of the zipped file
    """

    fileName = str(pid) + '_test_' + datetime.datetime.utcnow().strftime('%Y-%m-%d_%H-%M-%S') + '.zip'
    # Create the zip
    zip = zipfile.ZipFile(fileName, 'w')

    # Add individual files
    zip.write(str(pid) + "_dapt_param_settings.txt", compress_type=zipfile.ZIP_DEFLATED)

    # Add files in output
    for folder, subfolders, files in os.walk('output/'):
        for file in files:
            if ".png" not in file:
                zip.write(os.path.join(folder, file), 'output/'+file, compress_type = zipfile.ZIP_DEFLATED)

    # Add programming files
    zip.write('config/PhysiCell_settings.xml', compress_type = zipfile.ZIP_DEFLATED)
    zip.write('main*.cpp', compress_type = zipfile.ZIP_DEFLATED)
    zip.write('Makefile', compress_type = zipfile.ZIP_DEFLATED)
    zip.write('custom_modules/*', compress_type = zipfile.ZIP_DEFLATED)

    zip.close()

    return fileName