from __future__ import division
import numpy as np
import traceback as tb

dims = 3
cgRate = 8 
nParts = 40
qualNm = 'large' 
inName = 'N=' + str(nParts)
inPath = 'data/' + inName + '/' + qualNm + '/'
inName = inName + '_' + str(cgRate) + 'to1_' + qualNm
outNme = 'ref' + inName + '.trj'
outPat = 'data/'
frmNum = 100
trjNum = 10 

outp = open(outPat + outNme, "w")
outp.write(
"""1         999.000   999.000   999.000
1
0
4      3.0000
0
2      5.0000 0.9000
3
1.000
""")
outp.write('%5d    1\n' %(trjNum * frmNum))
outp.close()

for i in range(frmNum * trjNum):
    posLst = []
    velLst = []
    for k in range(dims):
        posLst.insert(k, [])
        velLst.insert(k, [])
    posTot = []
    velTot = []
    for k in range(dims):
        posTot.insert(k, 0.0)
        velTot.insert(k, 0.0)
    name = inPath + str(i + 1) + ".gro"
    print('Opening file ' + name)
    inpt = open(name, "r")
    for row in inpt:
        line = row.split()
        if(line == []):
            break
        elif(line[0] == 'Generated'):
            time = float(line[6])
            print(time)
        elif(len(line) == 1):
            atomNo = int(line[0])
            print(atomNo)
        elif(len(line) >= (2 + 2 * dims)):
            for k in range(dims):
                posLst[k].append(float(line[len(line) - 2 * dims + k]))
                velLst[k].append(float(line[len(line) - dims + k]))
        else:
            break
    if(len(posLst[0]) != atomNo):
        print('ERROR: there are %d data point' %len(posLst[0]))
    inpt.close()
    outp = open(outPat + outNme, "a")
    outp.write('%5.4f\n' %time)
    for j in range(atomNo):
        for k in range(dims):
            posTot[k] += posLst[k][j]
            velTot[k] += velLst[k][j]
        if(j % cgRate == (cgRate - 1)):
            for k in range(dims):
                outp.write('%5.4f ' %(posTot[k] / cgRate))
                posTot[k] = 0.0
            for k in range(dims):
                outp.write('%5.4f ' %(velTot[k] / cgRate))
                velTot[k] = 0.0
            outp.write('\n')
    outp.close()
