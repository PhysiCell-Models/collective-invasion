"""
    Example of how to use DAPT with Google Sheets.

    Before you can run this example, you need to create a new Google Sheet and add you Service Account Email as a collaborator with edit ability.  If you are still confused you can follow the documentation guide.
"""

import dapt
import os, sys, csv

sheet_id = '12dhlBi3BYBDa2oghU22ei0tUIbvtk5GS6T_sVOfify4'
db = dapt.sheets.Sheet(sheet_id, '../credentials.json')
ap = dapt.param.Param(db)

def init_sheet():
    header = ['id', 'startTime', 'endTime', 'status', 'a', 'b', 'c']
    for i in range(len(header)):
        db.sheet().update_cell(1, i+1, str(header[i]))

    data = [{'id':'t1', 'startTime':'2019-09-06 17:23', 'endTime':'2019-09-06 17:36', 'status':'finished', 'a':'2', 'b':'4', 'c':'6'}, {'id':'t2', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'10', 'c':''}, {'id':'t3', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'-10', 'c':''}]
    for i in range(3):
        for j in range(len(header)):
            db.sheet().update_cell(i+2, j+1, data[i][header[j]])

init_sheet()

input('If you go to the Google Sheet then you can see the sheet is now initialized.\nPress enter to run DAPT on the Google Sheet.')


while True:
    parameters = ap.next_parameters() #Get the next parameter
    if parameters == None:
        print("No more parameters to run!")
        break

    print("Request parameters: ")
    print(parameters)

    ap.successful(parameters["id"])


    try:
        ap.update_status(parameters['id'], 'Adding and inserting')

        c = int(parameters['a']) + int(parameters['b'])
        db.update_cell(db.get_row_index('id', parameters['id'])-2, 'c', c)

        # Update sheets to mark the test is finished
        ap.successful(parameters["id"]) #Test completed successfully so we mark it as such

        # End tests
    except ValueError:
        # Test failed
        print(ValueError)
        print("Test failed")
        ap.failed(parameters["id"], str(ValueError))

