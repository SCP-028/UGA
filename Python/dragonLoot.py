def addToInventory(inventory, addedItems):
    print('Inventory:')
    totalItem = 0
    for i in range(len(addedItems)):
        if addedItems[i] in inventory:
            inventory[str(addedItems[i])] += 1
        else:
            inventory[str(addedItems[i])] = 1
    for k, v in inventory.items():
        print(str(v) + ' ' + k)
        totalItem += int(v)
    print('Total number of items:' + str(totalItem))


inv = {'gold coin': 42, 'rope': 1}
dragonLoot = ['gold coin', 'dagger', 'gold coin', 'gold coin', 'ruby']
inv = addToInventory(inv, dragonLoot)
# if dragonLoot[0] in inv:
    # inv['gold coin'] += 1
# print(inv.items())
# print(len(dragonLoot))
