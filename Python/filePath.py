#! python3
"""
Chapter 8 - Reading and Writing Files.
"""
import os


def pathName():
    r"""
    Join path name. \ on Windows and / on Mac and Linux.
    """
    myFiles = ['accounts.txt', 'details.csv', 'invite.docx']
    for filename in myFiles:
        print(os.path.join('C:\\Users\\Joey', filename))


def currentDirectory():
    r"""
    cwd == Current Working Directory.
    chdir == Change Directory.
    """
    print(os.getcwd())
    os.chdir('C:\\Windows\\System32')


def thePaths():
    r"""
    absolute path: begin with root folder;
    relative path: relative to the program's current working directory;
    .: this directory;
    ..: home directory;
    """


def newFolder():
    r"""
    os.makedirs()
    """
