import random
def getAnswer(answerNumber):
        if answerNumber==1:
                return "It's certain"
        elif answerNumber==2:
                return "It's decidedly so"
        elif answerNumber==3:
                return "Yes"
        elif answerNumber==4:
        	return "Reply hazy try again"
        elif answerNumber==5:
        	return "Try again later"
r=random.randint(1,5)
fortune=getAnswer(r)
print(fortune)
