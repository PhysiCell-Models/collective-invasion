'''
    Ben Duggan
    1/18/19
    Create a CSV instance which can be used by param to get and run param sets
'''

import csv, os
class CSV:
    def __init__(self, csv_file, delimiter=','):
        self.csv_file = csv_file
        self.delimiter = delimiter

    def sheet(self):
        sheet = []
        with open(self.csv_file) as csvfile:
        	reader = csv.DictReader(csvfile)
        	for row in reader:
        		sheet.append(row)
        return sheet

    def getRecords(self):
        return self.sheet()

    def update_cell(self, i, j, text):
        i += 1 # Increment to account for header
        f = open(self.csv_file, 'r').readlines()

        assert i < len(f)
        
        # Find j index
        colNames = f[0][:-1].split(self.delimiter) #Get headers
        for r in range(len(colNames)):
            if j == colNames[r]:
                j = r
                break

        # Edit row
        row = f[i].split(self.delimiter)
        row[j] = text
        # Add new line if last header
        if j == len(colNames)-1:
            row[j] += '\n'

        # Edit CSV
        f[i] = ''
        for r in row:
            f[i] += r + self.delimiter
        f[i] = f[i][:-len(self.delimiter)]

        # Save CSV
        w = open(self.csv_file, 'w')
        for r in f:
            w.write(r)


if __name__ == '__main__':
    print('CSV sheet')
    csv = CSV('test.csv', ',')
    print(csv.sheet())
    print(csv.getRecords())
    csv.update_cell(0,1,'555')
    print(csv.getRecords())