'''
    Ben Duggan
    10/30/18
    Script for requesting and managing parameters
'''

import sheet
import datetime

class autoParam:
    def __init__(self):
        self.path = "DistributedAutomaticParameterTesting/"

        f = open(self.path+'config.txt', 'r').readlines()
        self.config = {}

        for i in range(0, len(f)):
            f[i] = f[i].replace("\n", "")
            self.config[f[i].split(":")[0]] = f[i].split(":")[1]
        print(self.path+"config.txt: ", self.config)

    def requestParameters(self):
        records = sheet.getRecords()

        if self.config["lastTest"] != "None":
            for i in range(0, len(records)):
                if int(self.config["lastTest"]) == records[i]["id"]:
                    records[i]["status"] = "in progress"
                    records[i]["start time"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                    records[i]["performed by"] = self.config["userName"]

                    sheet.sheet.update_cell(i+2, 2, "in progress")
                    sheet.sheet.update_cell(i+2, 3, records[i]["start time"])
                    sheet.sheet.update_cell(i+2, 5, self.config["userName"])

                    return records[i]

        for i in range(0, len(records)):
            if len(records[i]["status"]) == 0:
                records[i]["status"] = "in progress"
                records[i]["start time"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                records[i]["performed by"] = self.config["userName"]

                sheet.sheet.update_cell(i+2, 2, "in progress")
                sheet.sheet.update_cell(i+2, 3, records[i]["start time"])
                sheet.sheet.update_cell(i+2, 5, self.config["userName"])

                # Save id to local cache
                self.changeConfigLine("lastTest:None", "lastTest:"+str(records[i]["id"]))

                return records[i]

        return None

    def updateStatus(self, id, status):
        records = sheet.getRecords()
        index = -1

        # Remove id from local cache

        for i in range(0, len(records)):
            if records[i]["id"] == id:
                index = i

        if index == -1:
            return None

        records[index]["status"] = status

        sheet.sheet.update_cell(index+2, 2, status)

        return records[index]

    def parameterSuccessful(self, id):
        records = sheet.getRecords()
        index = -1

        # Remove id from local cache
        self.changeConfigLine("lastTest:"+str(id), "lastTest:None")

        for i in range(0, len(records)):
            if records[i]["id"] == id:
                index = i

        if index == -1:
            return None

        records[index]["status"] = "finished"
        records[index]["end time"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        sheet.sheet.update_cell(index+2, 2, "finished")
        sheet.sheet.update_cell(index+2, 4, records[index]["end time"])

        return records[index]

    def parameterFailed(self, id):
        records = sheet.getRecords()
        index = -1

        for i in range(0, len(records)):
            if records[i]["id"] == id:
                index = i

        if index == -1:
            return None

        records[index]["status"] = ""
        records[index]["start time"] = ""
        records[index]["performed by"] = self.config["userName"]
        records[index]["comments"] += " failed;"

        sheet.sheet.update_cell(index+2, 2, '')
        sheet.sheet.update_cell(index+2, 3, '')
        sheet.sheet.update_cell(index+2, 6, records[index]["comments"])

        return records[index]

    def checkForDBErrors(self):
        records = sheet.getRecords()

        for i in range(0, len(records)):
            # or records[i]["status"] == "run" or records[i]["status"] == "matlab" or records[i]["status"] == "ffmpeg" or records[i]["status"] == "zip" or records[i]["status"] == ""
            if (records[i]["status"] == "in progress") and len(records[i]["start time"]) > 0 and ((datetime.datetime.now()-datetime.datetime.strptime(records[i]["start time"], '%Y-%m-%d %H:%M:%S')).total_seconds() > int(self.config["resetTime"])):
                print("\"", records[i]["id"], "\" hasn't been marked as complete after running for: ", int((datetime.datetime.now()-datetime.datetime.strptime(records[i]["start time"], '%Y-%m-%d %H:%M:%S')).total_seconds()), " seconds. It has been marked as still needing to be ran.")

                records[i]["comments"] += "test possibly crashed"
                sheet.sheet.update_cell(i+2, 2, '')
                sheet.sheet.update_cell(i+2, 3, '')
                sheet.sheet.update_cell(i+2, 6, records[i]["comments"])

    def changeConfigLine(self, original, new):
        f = open(self.path+'config.txt', 'r').readlines()
        data = ""
        for i in range(0, len(f)):
            if f[i] == original+"\n":
                f[i] = new+"\n"
            data += f[i]

        with open(self.path+'config.txt', 'w') as file:
            file.writelines(data)
    def changeConfig(self):
        f = open(self.path+'config.txt', 'r').readlines()
        data = ""
        for i in self.config:
            data += i + ':' + self.config[i] + '\n'

        with open(self.path+'config.txt', 'w') as file:
            file.writelines(data)