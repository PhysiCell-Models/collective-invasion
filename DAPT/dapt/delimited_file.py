"""
    Delimited file
    ==============

    Create a CSV instance which can be used by param to get and run param sets
"""

import csv, os
from . import database

class Delimited_file(database.Database):
    """
        An interface for accessing and setting paramater set data.  
    """
        
    def __init__(self, csv_file, delimiter=','):
        
        super().__init__()
        
        self.csv_file = csv_file
        self.delimiter = delimiter

    def get_table(self):
        """
            Get the table from the database.

            Returns:
                An array with each element being a dictionary of the key-value pairs for the row in the database.
        """

        sheet = []
        with open(self.csv_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=self.delimiter)
            for row in reader:
                sheet.append(row)
        return sheet

    def get_keys(self):
        """
            Get the keys of the paramater set
            
            Returns:
                Array of strings with each element being a key (order is preserved if possible)
        """

        with open(self.csv_file, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=self.delimiter)
            return next(reader)

    def update_row(self, row_index, values):
        """
            Get the keys of the paramater set

            Args:
                row_index (int): the row id to replace
                values (OrderedDict): the key-value pairs that should be inserted
            
            Returns:
                A boolean that is True if successfully inserted and False otherwise.
        """

        header = self.get_keys()
        table = self.get_table()

        with open(self.csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=header)
            writer.writeheader()

            table[row_index] = values

            for row in table:
                writer.writerow(row)

            return True

    def update_cell(self, row_index, key, value):
        """
            Get the keys of the paramater set

            Args:
                row_index (int): the row id to replace
                key (str): the key of the value to replace
                value (str): the value to insert into the cell
            
            Returns:
                A boolean that is True if successfully inserted and False otherwise.
        """

        header = self.get_keys()
        table = self.get_table()

        with open(self.csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=header)
            writer.writeheader()

            table[row_index][key] = value

            for row in table:
                writer.writerow(row)

            return True 
            
    def get_row_index(self, column_key, row_value):
        """
            Get the row index given the column to look through and row value to match to.

            Args:
                column_key (str): the column to use.
                row_value (str): the row value to match with in the file and determin the row index.
            
            Returns:
                The index or -1 if it could not be determined
        """

        table = self.get_table()

        index = 0
        for row in table:
            if row[column_key] == row_value:
                return index
            index += 1

        return -1

if __name__ == '__main__':
    print('CSV sheet:')
    # Reset csv.  This is just used update the csv and reset it so interesting things happen.
    with open('test.csv', 'w') as f:
        writer = csv.DictWriter(f, fieldnames=['id', 'startTime', 'endTime', 'status', 'a', 'b', 'c'])
        writer.writeheader()
        writer.writerow({'id':'t1', 'startTime':'2019-09-06 17:23', 'endTime':'2019-09-06 17:36', 'status':'finished', 'a':'2', 'b':'4', 'c':'6'})
        writer.writerow({'id':'t2', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'10', 'c':''})
        writer.writerow({'id':'t3', 'startTime':'', 'endTime':'', 'status':'', 'a':'10', 'b':'-10', 'c':''})
    c = Delimited_file('test.csv', ',')
    print("get_keys: " + str(c.get_keys()))
    print("get_table: " + str(c.get_table()))
    print("update_cell: " + str(c.update_cell(1, "endTime", "09/02/19 12:10:00")))
    print("get_table: " + str(c.get_table()))
    print("update_row: " + str(c.update_row(1, {"id":"st-2", "startTime":"09/01/19 22:05:00", "endTime":""})))
    print("get_table: " + str(c.get_table()))