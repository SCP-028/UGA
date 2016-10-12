# ! python3
import sys
import pyperclip
"""
pw.py - An insecure password locker program.
"""
# Program Design and Data Structures
PASSWORDS = {'email': 'F7minlBDDuvMJuxESSKHFhTxFtjVB6',
             'blog': 'VmALvQyKAxiVH5G8v01if1MLZF3sdt',
             'luggage': '12345'}

# Handle Command Line Arguments
if len(sys.argv) < 2:
    print('Usage: py pw.py [account] - copy account password')
    sys.exit()
#  实际使用的时候输入 pw email就行，不用方括号
account = sys.argv[1]  # first command line arg is the account name

# Copy the Right Password
if account in PASSWORDS:
    pyperclip.copy(PASSWORDS[account])
    print('Password for ' + account + ' copied to clipboard.')
else:
    print('There is no account named ' + account)
