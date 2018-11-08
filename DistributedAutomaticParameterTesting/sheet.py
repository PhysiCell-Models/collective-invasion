# Ben Duggan
# 10/12/18
# Script for interfacing with Google Sheet

import time
import gspread
from oauth2client.service_account import ServiceAccountCredentials

spreedsheetID = "17QJFFXto0MbOX5dH9GFP3NevNHiuxj7eZf7Pevcg96U"

# use creds to create a client to interact with the Google Drive API
scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
creds = ServiceAccountCredentials.from_json_keyfile_name('credentials.json', scope)

def sheet():
    client = gspread.authorize(creds)

    return client.open_by_key(spreedsheetID).sheet1

def getRecords():
    return sheet().get_all_records()

if __name__ == '__main__':
    print(getRecords())

    time.sleep(15)

    print(getRecords())