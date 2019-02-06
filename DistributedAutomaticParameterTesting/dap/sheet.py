'''
    Ben Duggan
    1/8/19
    Main script to run distributed parameter testing
'''

import gspread
from oauth2client.service_account import ServiceAccountCredentials

class Sheet:
    def __init__(self, spreedsheetID, creds_path="credentials.json"):
        self.spreedsheetID = spreedsheetID
        self.scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
        self.creds = ServiceAccountCredentials.from_json_keyfile_name(creds_path, self.scope)

    # use creds to create a client to interact with the Google Drive API
    def sheet(self, sheet_id=0):
        client = gspread.authorize(self.creds)

        return client.open_by_key(self.spreedsheetID).get_worksheet(sheet_id)

    def getRecords(self):
        return self.sheet().get_all_records()

    def update_cell(self, i, j, text):
        if type(j) == type('a'):
            j = self.getKeyIndex(j)
        self.sheet().update_cell(i+2, j, text)

    def getKeyIndex(self, key):
        key_map = {}
        key_row = self.sheet().row_values(1)
        for i in range(len(key_row)):
            if key_row[i] == key:
                return i+1

if __name__ == '__main__':
    sheet = Sheet('1xZAbN6cs-89htm6EXkEYldQrSitzf5EnCGKwNl0a0Wo')
    print('Sheet result:')
    print(sheet.getRecords())
    print('Get key index:')
    print(sheet.getKeyIndex('id'))