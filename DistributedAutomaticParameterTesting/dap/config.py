'''
    Ben Duggan
    2/1/19
    Main script to run distributed parameter testing
'''

class Config:
	def __init__(self, path):
		self.path = path
		self.config = self.readConfig(self.path)

	def read_config(self):
		self.config = Config.readConfig(self.path)

	def change_config(self, key, value):
		self.config = Config.changeConfig(self.path, key, value)

	@staticmethod
	def readConfig(path):
	    f = open(path, 'r').readlines()
	    config = {}
	    for i in range(0, len(f)):
	        f[i] = f[i].replace("\n", "")
	        config[f[i].split(":")[0]] = f[i].split(":")[1]
	    return config

	'''
	@staticmethod
	def change_config(path, object):
		pass
	'''

	@staticmethod
	def changeConfig(path, key, value):
	    f = open(path, 'r').readlines()
	    data = ""
	    for i in range(0, len(f)):
	        if f[i].split(':')[0] == str(key):
	            f[i] = f[i].split(':')[0] + ':' + str(value) + '\n'
	        data += f[i]

	    with open(path, 'w') as file:
	        file.writelines(data)
	    return Config.readConfig(path)

class CommandLineArgs:
	@staticmethod
	def read():
		pass

if __name__ == '__main__':
	print(Config.readConfig('config.txt'))
	config = Config('config.txt')
	config.change_config('lastTest', 'k')
	print(config.config)