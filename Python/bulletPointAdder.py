# ! python3
"""
bulletPointAdder.py - Adds Wikipedia bulletpoints to the start
of each line of text on the clipboard.
"""
# Copy and Paste from the Clipboard
import pyperclip
text = pyperclip.paste()
pyperclip.copy(text)

# Separate the Lines of Text and Add the Star
lines = text.split('\n')
for i in range(len(lines)):
    lines[i] = '* ' + lines[i]

# Join the Modified Lines
text = '\n'.join(lines)
pyperclip.copy(text)
