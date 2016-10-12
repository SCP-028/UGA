#! python3
"""
tablePrinter.py - Takes a list of lists of strings and displays it
in a well-organized table with each column right-justified.
"""

def printTable(tableData):
#  TARGET RESULT
#  apples Alice  dogs
# oranges   Bob  cats
#cherries Carol moose
# banana David goose
    colWidths = [0] * len(tableData)
    for i in range(len(tableData)):
        colWidths[i] = len(max(tableData[i]))
    for j in range(len(tableData[0])):
       for k in range(len(tableData)):
        print(tableData[k][j].rjust(colWidths[k]),end = ' ')
        if k == len(tableData)-1:
            print()

tableData = [['apples','oranges','cherries','banana'],
            ['Alice','Bob','Carol','David'],
            ['dogs','cats','moose','goose'],
            ['asdfsadf','asdfsdgeg','wgsfdvbx','awrbfsdbdz']]
printTable(tableData)
