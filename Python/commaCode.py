spam = ['apples', 'bananas', 'tofu', 'cats']


def comma(list):
    for i in range(0, len(spam)-1):
        print(list[i], end=', ')
    print('and ' + list[len(list)-1])
comma(spam)
