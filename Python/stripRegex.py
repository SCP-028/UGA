#! python3
"""
Regex version of strip()
"""
import re
print('Please input the string to strip: ', end='')
stripString = input()
stripRegex = re.compile(r'\s*')
print(stripRegex.sub('', stripString))
