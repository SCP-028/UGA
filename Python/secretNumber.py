#This is a guess the number game.
import random
secretNumber = random.randint(1,20)
print('I am thinking of a number between 1 and 20.')

#Ask the player to guess 6 times.
for guessesTaken in range(1,7):
    print('Take a guess.')
    guess = int(input())
    if guess < secretNumber:
        print('Your guess is too low. ', end='')
    elif guess >secretNumber:
        print('Your guess is too high. ', end='')
    else:
        break     #This is the correct guess.
if guess == secretNumber:
    print('Good job! You got my number in ' + str(guessesTaken) + ' guesses!')
else:
    print('Nope. The number I was thinking was ' + str(secretNumber))
    
