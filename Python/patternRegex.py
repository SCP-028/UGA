#! python3
r"""
=====正则表达式=====
search(): 返回第一个匹配的字符串。
findall(): 返回所有匹配字符串。
    没有group ()返回list，有group返回tuple。
\d: 数字  \D: 非数字
\w: 字母、数字或_
\s: 空格、tab或回车
[]: 自定义表达式
{n}: 重复n次   {a,b}:重复a-b次    {n,}:n或更多
|: '或'
?: 可选模式，optional part of the pattern, 0 or 1.
*: 0 or n.
+: 1 or n.（至少出现一次！）
贪婪模式：默认，匹配最长字符串。非贪婪模式要在后面加?。
^: 在正则表达式开头，表明由后面的字符串开始。
$: 在正则表达式结尾，表明由后面的字符串结束。
.: 通配符（除回车外），匹配一个字符。 （re.DOTALL)
.*: 任意长度通配符。（贪婪模式）
re.I(GNORECASE): 不区分大小写
re.VERBOSE: 无视空格和注释
sub(): 取代字符串

"""
import re
def splitLine():
    """分割线"""
    print('==========')
def phoneNumber():
    r"""\d和{n}的使用"""
    phoneNumRegex = re.compile(r'(\d{3})-(\d{3}-\d{4})')
    mo = phoneNumRegex.search('My phone number is 400-882-3823')
    print('Phone number found: ' + mo.group())#group()返回整体结果，group(n)返回第n个括号结果
    print('Area code: ' + mo.group(1) + '\n' + 'Main number: ' + mo.group(2))
    splitLine()
def batMan():
    """|的使用"""
    batRegex = re.compile(r'Bat(man|mobile|copter|bat)')
    heroRegex = re.compile(r'Batman|Tina Fey')
    mo = batRegex.search('Batmobile lost a wheel.')
    mo1 = heroRegex.search('Batman and Tina Fey')
    mo2 = heroRegex.search('Tina Fey and Batman')#两个目标都存在，第一个作为结果。
    print(mo.group())
    print(mo.group(1))
    print(mo1.group())
    print(mo2.group())
    splitLine()
def batWoman():
    """?的使用"""
    batRegex = re.compile(r'Bat(wo)?man')
    mo1 = batRegex.search('The Adventures of Batman')
    mo2 = batRegex.search('The Adventures of Batwoman')
    print(mo1.group())
    print(mo2.group())
    splitLine()
def batWowoman():
    """
    *的使用，+类似。
    """
    batRegex = re.compile(r'Bat(wo)*man')
    mo1 = batRegex.search('The Adventures of Batman')
    mo2 = batRegex.search('The Adventures of Batwoman')
    mo3 = batRegex.search('The Adventures of Batwowowowowoman')
    print(mo1.group() + '\n' + mo2.group() + '\n' + mo3.group())
    splitLine()
def greedyMode():
    """贪婪模式和非贪婪模式"""
    greedyHaRegex = re.compile(r'(Ha){3,5}')
    nongreedyHaRegex = re.compile(r'(Ha){3,5}?')
    mo1 = greedyHaRegex.search('HaHaHaHaHa')
    mo2 = nongreedyHaRegex.search('HaHaHaHaHa')
    print(mo1.group() + '\n' + mo2.group())
    splitLine()
def xMas():
    r"""findall(),\s,\w"""
    xmasRegex = re.compile(r'\d+\s\w+')
    print(xmasRegex.findall('12 drummers, 11 pipers, 10 lords, 9 ladies, '+\
        '8 maids, 7 swans, 6 geese, 5 rings, 4 birds, 3 hens, '+\
        '2 doves, 1 partridge'))
    splitLine()
def findVowel():
    r"""[]内为要的；[^xxx]是反向匹配（不要xxx）;
    [a-zA-Z0-9.]匹配所有字母和数字和.(不用\.)
    """
    vowelRegex = re.compile(r'[aeiouAEIOU]')
    print(vowelRegex.findall('Robocop eats baby food. BABY FOOD.'))
    splitLine()
def subRegex():
    """sub()的使用"""
    namesRegex = re.compile(r'Agent \w+')
    print(namesRegex.sub('UNKNOWN', 'Agent Alice gave the secret documents to Agent Bob.'))
    agentRegex = re.compile(r'Agent (\w)\w*')
    print(agentRegex.sub(r'Agent \1****', 'Agent Alice told Agent Carol'\
     ' that Agent Eve knew Agent Bob was a double agent.'))
    splitLine()
def phoneRegex():
    """re.VERBOSE的使用"""
    phonesRegex = re.compile(r'''(
        (\d{3}|\(\d{3}\))?      #area code
        (\s|-|\.)?          #separator
        \d{3}               #first 3 digits
        (\s|-|\.)           #separator
        \d{4}               #last 4 digits
        (\s*(ext|x|ext.)\s*\d{2,5})?   #extension
        )''', re.VERBOSE)
    print(phonesRegex.search('The phone number is 400-882-3823').group())
phoneNumber()
batMan()
batWoman()
batWowoman()
greedyMode()
xMas()
findVowel()
subRegex()
phoneRegex()
