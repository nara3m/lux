#!/usr/bin/env python
# Usage : <program_name> <input_file> <output_file> <log_file>

from subprocess import Popen, PIPE
import sys, re, string, warnings, time

pipe = Popen("date", stdout=PIPE)
start_date = pipe.communicate()[0]
ticks = time.time()

output_file = open(sys.argv[2], 'w')

# sub-routine for calculating levenshtein distance of two strings 
# taken from http://code.activestate.com/recipes/576874-levenshtein-distance/

#def printMatrix(m):
    #print ' '
    #for line in m:
        #spTupel = ()
        #breite = len(line)
        #for column in line:
            #spTupel = spTupel + (column, )
        #print "%3i"*breite % spTupel

def levenshtein(each_SMILES_string, each_SMILES_string_again):
	l1 = len(each_SMILES_string)
	l2 = len(each_SMILES_string_again)

	matrix = [range(l1 + 1)] * (l2 + 1)
	for zz in range(l2 + 1):
		matrix[zz] = range(zz,zz + l1 + 1)
	for zz in range(0,l2):
		for sz in range(0,l1):
			if each_SMILES_string[sz] == each_SMILES_string_again[zz]:
				matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz])
			else:
				matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz] + 1)
			#print "That's the Levenshtein-Matrix:"
			#printMatrix(matrix)
	return matrix[l2][l1]
# end of sub-routine

iterator2=0

i = open(sys.argv[1], 'r')
for each_line in i:

	iterator1=0
	split_each_line_into_two = string.split(each_line)
	output_file.write(split_each_line_into_two[1]+",")
	
	j = open(sys.argv[1], 'r')
	for each_line_again in j:
		
		if iterator1==iterator2:

			output_file.write("0\n")				
			break

		else:
			split_each_line_into_two_again = string.split(each_line_again)

			distance = levenshtein(split_each_line_into_two[0], split_each_line_into_two_again[0])
			max_length = float(max(len(split_each_line_into_two[0]), len(split_each_line_into_two_again[0])))
			distance = distance/max_length
			distance = "{0:.4f}".format(distance)
			output_file.write(str(distance)+",")
		
		iterator1=iterator1+1
			
	j.close()
	iterator2=iterator2+1
	
i.close()
output_file.close()

pipe = Popen("date", stdout=PIPE)
end_date = pipe.communicate()[0]
ticks2 = time.time()-ticks

# log file
# f4=open(sys.argv[3],'a')
# f4.write("\nMethod : Levenshtein Distance\nInput file : "+sys.argv[1]+"\nOutput file : "+sys.argv[2]+"\nStart time : "+start_date+"End time : "+end_date+"Time elapsed : "+str(ticks2)+" seconds\n")
# f4.close()
