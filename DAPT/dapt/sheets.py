"""
Sheets
====== 

Class which allows for Google Sheets to be used as paramater set database.



Note: if you have data in the first row, you must have entries in some other row.
"""

import gspread
from oauth2client.service_account import ServiceAccountCredentials
from . import database

class Sheet(database.Database):
    """
    An interface for accessing and setting paramater set data.  You must either provide a Config object or client_id and client_secret.

    Keyword Args:
        config (Config): A Config object which contains the client_id and client_secret. 
        spreedsheet_id (str): the Google Sheets ID
        creds (str): the path to the file containing the Google API credentials.  Default is ``credentials.json``.
        sheet_id (int): the the sheet id to use.  0 is used if no value is givin for sheet_title, sheet_id or in the Config
        sheet_title (str): the title of the sheet to use
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__()
        
        self.spreedsheetID = None
        self.scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
        self.creds = None
        self.sheet_id = -1
        self.sheet_title = None
        self.config = None

        if len(kwargs) == 0:
            raise ValueError("You must provide a Config object or spreedsheetID and credentials file.")
        if 'config' in kwargs:
            self.config=kwargs['config']

            if 'sheets-spreedsheet-id' in self.config.config:
                self.spreedsheetID = self.config.config['sheets-spreedsheet-id']
            if 'sheets-creds-path' in self.config.config:
                self.creds = ServiceAccountCredentials.from_json_keyfile_name(self.config.config['sheets-creds-path'], self.scope)
            if 'sheets-worksheet-id' in self.config.config:
                self.sheet_id = self.config.config['sheets-worksheet-id']
            if 'sheets-worksheet-title' in self.config.config:
                self.sheet_title = self.config.config['sheets-worksheet-title']
        if not self.spreedsheetID:
            if 'spreedsheet_id' in kwargs:
                self.spreedsheetID = kwargs['spreedsheet_id']
            else:
                raise ValueError("You must specify the spreedsheet id in the arguments or config.")
        if not self.creds:
            if 'creds' in kwargs:
                self.creds = ServiceAccountCredentials.from_json_keyfile_name(kwargs['creds'], self.scope)
            else:
                self.creds = ServiceAccountCredentials.from_json_keyfile_name('credentials.json', self.scope)
        if not self.sheet_title and self.sheet_id == -1:
            if 'sheet_title' in kwargs:
                self.sheet_title = kwargs['sheet_title']
        if self.sheet_id == -1 and not self.sheet_title:
            if 'sheet_id' in kwargs:
                self.sheet_id = kwargs['sheet_id']
            else:
                self.sheet_id = 0

        self.client = gspread.authorize(self.creds)
        self.sheet = self.client.open_by_key(self.spreedsheetID)

    def _update_creds(self):
        """
        Update the credentials if they are expired
        """

        if self.creds.access_token_expired:
            self.client.login()

    def worksheet(self, *args, **kwargs):
        """
        Get a Google Sheet object.  The worksheet id or title are obtained from the Config file or initialization.

        Returns:
            A Google Sheet worksheet
        """

        self._update_creds()

        if self.sheet_title:
            print('sup dogg')
            return self.sheet.worksheet(self.sheet_title)
        elif self.sheet_id >= 0:
            return self.sheet.get_worksheet(self.sheet_id)
        else:
            return self.sheet.get_worksheet(0)

    def get_table(self):
        """
        Get the table from the database.

        Returns:
            An array with each element being a dictionary of the key-value pairs for the row in the database.
        """

        self._update_creds()

        return self.worksheet().get_all_records()

    def get_keys(self):
        """
        Get the keys of the paramater set
        
        Returns:
            Array of strings with each element being a key (order is preserved if possible)
        """

        self._update_creds()

        return self.worksheet().row_values(1)

    def update_row(self, row_index, values):
        """
        Get the row of the paramater set

        Args:
            row_index (int): the index of the row to replace (starting from 1).  Indices less than 1 will return False.  Indices greater than the table length will be appended.
            values (Dict): the key-value pairs that should be inserted.  If the dictionary contains more values then number of columns, the table will be extended.
        
        Returns:
            A boolean that is Trues if successfully inserted and False otherwise.
        """

        if row_index < 0:
            return False

        row = [[]]
        for i in values:
            row[0].append(values[i])
        
        self._update_creds()

        start = gspread.utils.rowcol_to_a1(row_index+2, 1)
        end = gspread.utils.rowcol_to_a1(row_index+2, len(values))
        range_label = '%s!%s:%s' % (self.worksheet().title, start, end)
        
        try:
            return self.sheet.values_update(range_label, params={'valueInputOption': 'RAW'}, body={'values': row})
        except:
            return False
        return True

    def update_cell(self, row_id, key, value):
        """
        Get the keys of the paramater set
        
        Args:
            row_id (int): the row id to replace
            key (str): the key of the value to replace
            value (str): the value to insert into the cell
        
        Returns:
            A boolean that is True if successfully inserted and False otherwise.
        """
        try:
            self._update_creds()
            self.worksheet().update_cell(row_id+2, self.get_key_index(key)+1, str(value))
        except Exception as e:
            raise e
        return True

    def get_key_index(self, column_key):
        """
        Get the column index given the key.

        Args:
            column_key (str): the key to find the index of
        
        Returns:
            The index or -1 if it could not be determined.
        """

        self._update_creds()

        key_map = {}
        key_row = self.worksheet().row_values(1)
        for i in range(len(key_row)):
            if str(key_row[i]) == str(column_key):
                return i
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

        self._update_creds()

        col = self.worksheet().col_values(self.get_key_index(column_key)+1)
        for i in range(len(col)):
            if str(col[i]) == str(row_value):
                return i-1
        return -1
    