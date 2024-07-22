def po(list0):
    list1 = list0
    for i in range(len(list0)):
        if list0[i] < 0:
            list1[i] = 0
    return list1