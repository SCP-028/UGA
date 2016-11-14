from collections import Counter
inventory = Counter({'gold': 42,
                    'rope': 1})


def addToInventory(inventory, addedItems):
    inventory.update(addedItems)


addedItems = ['gold', 'gold', 'dagger', 'ruby']
addToInventory(inventory, addedItems)
totalNum = 0
for k, v in inventory.items():
    print(k + ': ' + str(v))
    totalNum += int(v)
print('Total items: ' + str(totalNum))
