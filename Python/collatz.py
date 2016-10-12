def collatz(number):
    print(number)
    while number != 1:
        if number % 2 == 0:
            number = number//2
            print(number)
        elif number % 2 == 1:
            number = 3*number+1
            print(number)
print('Please type a number:', end='')
try:
    number = int(input())
    collatz(number)
