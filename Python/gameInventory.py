things = {'arrow': 12, 'gold coin': 42, 'rope': 1, 'torch': 6, 'dagger': 1}


def displayInventory(inventory):
    print('Inventory: ')
    numThings = 0
    for k,  v in inventory.items():
        print(str(v) + ' ' + k)
        numThings += int(v)
    print('Total number of items:  ' + str(numThings))
displayInventory(things)
