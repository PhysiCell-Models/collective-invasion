'''
    Ben Duggan
    1/18/19
    Main script to run distributed parameter testing
'''

import datetime
#import config

class Param:
    def __init__(self, config, sheet):
        self.path = "DistributedAutomaticParameterTesting/"

        self.sheet = sheet

        self.config = config
        self.conf = config.config

        if self.conf['numOfRuns']:
            self.count = 0
        self.lastParam = None

    def requestParameters(self):
        if self.conf['numOfRuns']:
            if int(self.conf['numOfRuns']) == -1 or self.count < int(self.conf['numOfRuns']):
                    self.count += 1
            else:
                return None

        records = self.sheet.getRecords()

        if "lastTest" in self.conf and self.conf["lastTest"] != "None":
            print("Using lastTest from config.txt")
            for i in range(0, len(records)):
                if str(self.conf["lastTest"]) == str(records[i]["id"]):
                    if "startTime" in records[i]:
                        records[i]["startTime"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                        self.sheet.update_cell(i, 'startTime', records[i]["startTime"])
                    if "performedBy" in records[i]:
                        records[i]["performedBy"] = self.conf["userName"]
                        self.sheet.update_cell(i, 'performedBy', self.conf["userName"])

                    self.lastParam = records[i]["id"]

                    return records[i]

        for i in range(0, len(records)):
            if len(records[i]["status"]) == 0 and (self.lastParam == None or (not str(records[i]["id"]) == str(self.lastParam))):
                try:
                    if int(self.conf['computerStrength']) < int(records[i]["computerStrength"]):
                        continue
                except:
                    pass

                if "status" in records[i]:
                    records[i]["status"] = "in progress"
                    self.sheet.update_cell(i, 'status', "in pogress")
                if "startTime" in records[i]:
                    records[i]["startTime"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                    self.sheet.update_cell(i, 'startTime', records[i]["startTime"])
                if "performedBy" in records[i]:
                    records[i]["performedBy"] = self.conf["userName"]
                    self.sheet.update_cell(i, 'performedBy', self.conf["userName"])

                self.lastParam = records[i]["id"]

                # Save id to local cache
                self.config.change_config("lastTest", str(records[i]["id"]))

                return records[i]

        return None

    def updateStatus(self, id, status):
        records = self.sheet.getRecords()
        index = -1

        # Remove id from local cache

        for i in range(0, len(records)):
            if str(records[i]["id"]) == str(id):
                index = i

        if index == -1:
            return None

        records[index]["status"] = status
        self.sheet.update_cell(index, 'status', status)

        return records[index]

    def parameterSuccessful(self, id):
        records = self.sheet.getRecords()
        index = -1

        # Remove id from local cache
        self.config.change_config("lastTest", "None")

        for i in range(0, len(records)):
            if str(records[i]["id"]) == str(id):
                index = i

        if index == -1:
            return None

        if 'status' in records[index]:
            records[index]["status"] = "finished"
            self.sheet.update_cell(index, 'status', "finished")
        if 'endTime' in records[index]:
            records[index]["endTime"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            self.sheet.update_cell(index, 'endTime', records[index]["endTime"])

        return records[index]

    def parameterFailed(self, id):
        records = self.sheet.getRecords()
        index = -1

        for i in range(0, len(records)):
            if str(records[i]["id"]) == str(id):
                index = i

        if index == -1:
            return None

        if 'status' in records[index]:
            records[index]["status"] = ""
            self.sheet.update_cell(index, 'status', '')
        if 'endTime' in records[index]:
            records[index]["endTime"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            self.sheet.update_cell(index, 'endTime', records[index]['endTime'])
        if 'comments' in records[index]:
            records[index]["comments"] += " failed;"
            self.sheet.update_cell(index, 'comments', records[index]["comments"])

        return records[index]

    def checkForDBErrors(self):
        records = self.sheet.getRecords()

        if len(records) > 0:
            if 'startTime' not in records[1] or 'resetTime' not in self.conf:
                return None

        for i in range(0, len(records)):
            if len(records[i]["startTime"]) > 0 and ((datetime.datetime.now()-datetime.datetime.strptime(records[i]["startTime"], '%Y-%m-%d %H:%M:%S')).total_seconds() > int(self.conf["resetTime"])):
                print("\"", records[i], "\" hasn't been marked as complete after running for: ", int((datetime.datetime.now()-datetime.datetime.strptime(records[i]["startTime"], '%Y-%m-%d %H:%M:%S')).total_seconds()), " seconds. It has been marked as still needing to be ran.")

                if 'comments' in records[i]:
                    records[i]["comments"] += "test possibly crashed"
                    self.sheet.update_cell(i, 'comments', records[i]["comments"])
                if 'status' in records[i]:
                    self.sheet.update_cell(i, 'status', '')
                if 'startTime' in records[i]:
                    self.sheet.update_cell(i, 'startTime', '')
                return records
        return None    

