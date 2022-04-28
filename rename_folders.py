import json
import os

originalsPath = "input/"

projectFolder = ""

path = originalsPath + projectFolder

with open('samples.json', 'r') as f:
  samples = json.load(f)

for i in samples["samples"]:
	sampleName = samples["samples"][i]["name"]
	os.rename(path+i, path+sampleName)
	print("Renamed "+ i +" to " + sampleName)
