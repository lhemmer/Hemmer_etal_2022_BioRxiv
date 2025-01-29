import sys

fileName = sys.argv[1]
file = open(fileName, "r")
file = file.readlines()

newFileName = sys.argv[2]
newFile = open(newFileName, "w+")

dnFileName = sys.argv[3]
newFileDN = open(dnFileName, "w+")

libName = sys.argv[4]
lib = open(libName, "r")


#create dictionary to store contig library entries, descriptor + contig = key
libDict = {}

#enter contig library entries
for line in lib:
	line = line.strip()
	#skip header lines
	if(line[0] == '#' or line[0] == '"'):
		continue

	else:
		libLine = line.split(",")
		#fuse contig name and coordinates to handle contigs with heterogenous chromatin states (multiple entries)
		key = (libLine[3] + ";" + libLine[4])
		
		##for now, just make dict with contig, don't worry about contigs with both het and euchromatin
		#key = (libLine[3])

		#value = remaining info about each contig
		value = [libLine[0], libLine[1], libLine[2], libLine[4], libLine[5]]

		libDict[key] = value
		#print(libDict)


#write header line
newFile.write("contig" + "\t" + "chromosome" + "\t" + "start" + "\t" + "end" + "\t" + "ID" + "\t" + "called.by" + "\t" + "chromatin.state" + "\t" + "region" + "\t" + "ref" + "\t" + "population" + "\n")
newFileDN.write("contig" + "\t" + "chromosome" + "\t" + "start" + "\t" + "end" + "\t" + "ID" + "\t" + "called.by" + "\t" + "chromatin.state" + "\t" + "region" + "\t" + "ref" + "\t" + "population" + "\n")

#process output file from McClintock
for line in file:
	
	line = line.strip()
	fields = line.split("\t")
	
	#skip header line
	if(fields[0][0:5] == "track"):
		population = (line.split("=")[1]).split("_")[0]
		population = population[1:]
		continue
	

	else:
		#handle contigs with different regions of chromatin
		fileContig = fields[0]
		coordStart = int(fields[1])
		coordEnd = int(fields[2])

		#newFile.write("file contig: " + fileContig + "\n")
		libKeys = libDict.keys()
		
		for key in libKeys:
			keyFields = key.split(';')
			#newFile.write(keyFields[0] + " " + keyFields[1] + "\n")
			if(fileContig == keyFields[0]):
				#newFile.write("matched contig " + keyFields[0] + "\n")
				if(keyFields[1] == ""):
					#newFile.write(key + " is a homogenous contig" + "\n")
					actualKey = key	
					
				if(keyFields[1] != ""):
					keyStart = int(keyFields[1].split("-")[0])
					keyEnd = int(keyFields[1].split("-")[1])
					#newFile.write("keyStart: " + str(keyStart) + "\n")
					#newFile.write("keyEnd: " + str(keyEnd) + "\n")
					#newFile.write("coordStart: " + str(coordStart) + "\n")
					#newFile.write("coordEnd: " + str(coordEnd) + "\n")

					#this does not make sense
					#newFile.write("coordStart " + str(coordStart) + " > " + "keyStart " + str(keyStart) + " " + str(coordStart > keyStart) + "\n") 
					#newFile.write("coordEnd " + str(coordEnd) + " < " + "keyEnd " + str(keyEnd) + " " + str(coordEnd < keyEnd) + "\n")
					
					#newFile.write("debug" + "\n")
					#newFile.write(str(keyStart) + "\n")
					#newFile.write(str(keyEnd) + "\n")
					#newFile.write(str(coordEnd) + "\n")
					#newFile.write(str(coordStart) + "\n")

					if(keyStart < coordStart and keyEnd > coordEnd):
							#newFile.write("logic test")
							#newFile.write(key + " is a heterogenous contig" + "\n")
							actualKey = key
							#newFile.write("actual Key: " + actualKey + "\n")
					else:
						#newFile.write("else" + "\n")
						#newFile.write(key + " did not match " + fileContig + " " + str(coordStart) + "-" + str(coordEnd) + "\n")
						continue
						

		#get info about McC call by finding value for which contig is contained in key, and range matches too
		info = libDict[actualKey]
		
		#record which program called this TE
		#IDfields = fields[3].split("_")
		IDfields = fields[3]
		program = ""
		if("temp" in IDfields):
			program = "temp"
		
		if("ngs_te_mapper" in IDfields):
			program = "ngs_te_mapper"

		if("telocate" in IDfields):
			program = "telocate"

		if("popoolationte" in IDfields):
			program = "popoolationte"
		
		if("retroseq" in IDfields):
			program = "retroseq"
			
		if("popoolationte2" in IDfields):
			program = "popoolationte2"
		
		if("relocate2" in IDfields):
			program = "relocate2"
	

		#reference or non-reference
		if("non-reference" in fields[3]):
			#create string to be added to and written to output file
			#format: contig, pos, TE, program, chromatin state, bin/region
			deNovo = "non-reference"
			ID = fields[3].split("_non-reference")[0]
			#print processed output to new file
			newFileDN.write(fields[0] + "\t" + info[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + ID + "\t" + program + "\t" + info[1] + "\t" + info[2] + "\t" + deNovo + "\t" + population + "\n")
		
		else:
			#create string to be added to and written to output file
			#format: contig, pos, TE, program, chromatin state, bin/region
			deNovo = "reference"
			ID = fields[3].split("_reference")[0]

			if(info[2] == ""):
				info[2] = "."

			#print processed output to new file
			newFile.write(fields[0] + "\t" + info[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + ID + "\t" + program + "\t" + info[1] + "\t" + info[2] + "\t" + deNovo + "\t" + population + "\n")


				


		
		
