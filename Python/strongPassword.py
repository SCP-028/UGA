#! python3
"""
strongPassword.py
检测密码强度，必须8字符以上，有大写&小写字母，至少1个数字。
"""
import re
pass1Regex = re.compile(r'.{8}')
pass2Regex = re.compile(r'[A-Z]+')
pass3Regex = re.compile(r'[a-z]+')
pass4Regex = re.compile(r'\d+')


def passwordCheck():
    """Check the password."""
    print('Please input your password: ', end='')
    password = input()
    i = 0
    if pass1Regex.search(password):
        i += 1
    if pass2Regex.search(password):
        i += 2
    if pass3Regex.search(password):
        i += 3
    if pass4Regex.search(password):
        i += 4
    if i == 10:
        print('Nice password!')
    else:
        print('Check your pass!')
passwordCheck()
