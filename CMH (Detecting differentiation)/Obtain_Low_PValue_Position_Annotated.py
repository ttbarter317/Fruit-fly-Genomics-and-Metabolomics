inputFile = "meta_20cov_sigSNP.txt"
outputFileName = "outputmeta_20cov_sigSNP.txt"

#Section 1: Windows
#This sets up the windows in which you want to find the lowest p-value therin
window = {'X': [], '2L': [], '2R': [], '3L': [], '3R': []}
#This is probably the only number you need to modify in section 1
kbSize = 50000
num_windows = int (33000000 / kbSize)

for i in range(num_windows):
    window['X'].append([])
    window['2L'].append([])
    window['2R'].append([])
    window['3L'].append([])
    window['3R'].append([])
#End Section 1

#Section 2: Read File
#This section read the file of SNP positions and their p-values
with open(inputFile, 'r') as file:
    for line in file:
        line = line.split()
        #Depending on how your files is written up, you might need to change any of these numbers
        _chr = line[0]
        pos = (int)(line[1])
        #Except this one is one you will almost definitely need to change
        #Change to whatever column is the p-value column, but remember 0 based indexing is in effect
        p_value = (float)(line[64])
        out = [_chr, pos, p_value]
        window[_chr][(int)(pos / kbSize)].append(out)

#Section 3: The Test
#This section was there for the specific test condition for the FLAME paper
#Kept here for posterity
'''
extraInput = "output50kb.txt"
extraTest = {'X': [], '2L': [], '2R': [], '3L': [], '3R': []}
with open(extraInput, 'r') as extra:
    for line in extra:
        line = line.split()
        p_value = (float)(line[2])
        if p_value != -1:
            num = (int)(line[1])
            num = (int) (num / 50000)
            num = num * 50000
            for i in range(5):
                extraTest[line[0]].append(num + i * 10000)
'''

#Section 4: Finding low p-values and output
output = ""
for i in list(window.keys()):
    for j in range(num_windows):
        min_p = 9999
        best = [i, j * kbSize, -1]
        #This is probably the main part you will want to change
        #The test condition for the FLAME paper was there had to be at least
        #3 minimum positions inside a 50kb window before it would find the
        #lowest p-value in the 10kb window
        #if len(window[i][j]) > 2 or (j * kbSize) in extraTest[i]
        test = True
        if test:
            for k in window[i][j]:
                if k[2] < min_p:
                    min_p = k[2]
                    best = k
        #This if also probably unnecessary for 
        if best[2] != -1:
            output += best[0] + '\t' + str(best[1]) + '\t' + str(best[2]) + '\n'

with open(outputFileName, 'w') as outputFile:
    outputFile.write(output)
