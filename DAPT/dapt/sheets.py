"""
Sheets
====== 

Class which allows for Google Sheets to be used as paramater set database.
"""


import gspread
from oauth2client.service_account import ServiceAccountCredentials
from . import database

class Sheet(database.Database):
    """
        An interface for accessing and setting paramater set data.

        Args:
            spreedsheetID (str): the Google Sheets ID
            creds_file (str): the path to the file containing the Google API credentials
            sheet_id (int): the the sheet id to use (0 by default)
    """

    def __init__(self, spreedsheetID, creds_file="credentials.json", sheet_id=0):
        
        super().__init__()

        self.spreedsheetID = spreedsheetID
        self.scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
        self.creds = ServiceAccountCredentials.from_json_keyfile_name(creds_file, self.scope)
        self.sheet_id = sheet_id
    
    def sheet(self):
        """
            Get a Google Sheet object

            Returns:
                A Google Sheet worksheet
        """

        client = gspread.authorize(self.creds)

        return client.open_by_key(self.spreedsheetID).get_worksheet(self.sheet_id)

    def get_table(self):
        """
            Get the table from the database.

            Returns:
                An array with each element being a dictionary of the key-value pairs for the row in the database.
        """

        return self.sheet().get_all_records()

    def get_keys(self):
        """
            Get the keys of the paramater set
            
            Returns:
                Array of strings with each element being a key (order is preserved if possible)
        """

        return self.sheet().row_values(1)

    def update_row(self, row_id, values):
        """
            Get the row of the paramater set

            Args:
                row_id (int): the row id to replace
                values (OrderedDict): the key-value pairs that should be inserted
            
            Returns:
                A boolean that is True if successfully inserted and False otherwise.
        """

        for i in values:
            try:
                self.sheet().update_cell(row_id+2, self.get_key_index(i), str(values[i]))
            except:
                return False
        return True

    def update_cell(self, row_id, key, value):
        """
            Get the keys of the paramater set
            
            Args:
                row_id (str): the row id to replace
                key (str): the key of the value to replace
                value (str): the value to insert into the cell
            
            Returns:
                A boolean that is True if successfully inserted and False otherwise.
        """
        try:
            self.sheet().update_cell(row_id+2, self.get_key_index(key), str(value))
        except:
            return False
        return True

    def get_key_index(self, column_key):
        """
            Get the column index given the key.

            Args:
                column_key (str): the key to find the index of
            
            Returns:
                The index or -1 if it could not be determined.
        """

        key_map = {}
        key_row = self.sheet().row_values(1)
        for i in range(len(key_row)):
            if key_row[i] == column_key:
                return i+1
        return -1

    def get_row_index(self, column_key, row_value):
        """
            Get the row index given the column to look through and row value to match to.

            Args:
                column_key (str): the key to find the index of
                row_value (str): the value of the cell to fine
            
            Returns:
                The index or -1 if it could not be determined.
        """

        col = self.sheet().col_values(self.get_key_index(column_key))
        for i in range(len(col)):
            if col[i] == row_value:
                return i+1
        return -1
        

if __name__ == '__main__':
    sheet = Sheet('1xZAbN6cs-89htm6EXkEYldQrSitzf5EnCGKwNl0a0Wo')
    print('Sheet result:')
    print(sheet.getRecords())
    print('Get key index:')
    print(sheet.getKeyIndex('id'))