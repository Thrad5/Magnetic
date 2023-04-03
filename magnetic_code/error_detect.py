"""
Created on Wed Mar 1 15:20:32 2023

@author: ram86
"""

def inpTF(msg):
    while True:
        test = input(msg)
        if test == 'y':
            return True
        elif test == 'n':
            return False
        else:
            print('That is an invalid input please enter either y or n')

def inpInt(msg,minm = None, maxm = None):
    while True:
        if minm != '':
            if maxm != '':
                try:
                    test = int(input(msg))
                except ValueError:
                    print('The number must be an integer.')
                    continue
                if test < minm or test > maxm:
                    print(f'The number must be between {minm} and {maxm}.')
                else:
                    return test
            else:
                try:
                    test = int(input(msg))
                except ValueError:
                    print('The number must be an integer.')
                    continue
                if test < minm:
                    print(f'The number must be larger than {minm}.')
                else:
                    return test
        else:
            if maxm != '':
                try:
                    test = int(input(msg))
                except ValueError:
                    print('The number must be an integer.')
                    continue
                if test > maxm:
                    print(f'The number must be smaller than {maxm}.')
                else:
                    return test
            else:
                try:
                    test = int(input(msg))
                except ValueError:
                    print('The number must be an integer.')
                    continue
                else:
                    return test

def inpFlt(msg,minm = None, maxm = None):
    while True:
        if minm != None:
            if maxm != None:
                try:
                    test = float(input(msg))
                except ValueError:
                    print('The number must be an decimal.')
                    continue
                if test < minm or test > maxm:
                    print(f'The number must be between {minm} and {maxm}.')
                else:
                    return test
            else:
                try:
                    test = float(input(msg))
                except ValueError:
                    print('The number must be an decimal.')
                    continue
                if test < minm:
                    print(f'The number must be larger than {minm}.')
                else:
                    return test
        else:
            if maxm != None:
                try:
                    test = float(input(msg))
                except ValueError:
                    print('The number must be an decimal.')
                    continue
                if test > maxm:
                    print(f'The number must be smaller than {maxm}.')
                else:
                    return test
            else:
                try:
                    test = float(input(msg))
                except ValueError:
                    print('The number must be an decimal.')
                    continue
                else:
                    return test

def inpDate():
    while True:
        try:
            yr = int(input('What year do you want the magnetic field for? (1900-2030):'))
        except ValueError:
            print('The year must be an integer.')
            continue
        if yr < 1900 or yr > 2030:
            print('The year must be between 1900 and 2030.')
        else:
            break
    if yr % 400 == 0:
        leap = True
    elif yr % 100 == 0:
        leap = False
    elif yr % 4 == 0:
        leap = True
    else:
        leap = False
    while True:
        try:
            mth = int(input('What month do you want the magnetic field for as a number?'))
        except ValueError:
            print("The month must be a number.")
            continue
        if mth < 1 or mth > 12:
            print("There are only 12 months in a year.")
        else:
            break
    day30 = [4,6,9,11]
    day31 = [1,3,5,7,8,10,12]
    if mth in day30:
        day = inpInt('What day of the month do you want the magnetic field for? (1-30)',1,30)
    elif mth in day31:
        day = inpInt('What day of the month do you want the magnetic field for? (1-31)',1,31)
    elif leap == True:
        day = inpInt('What day of the month do you want the magnetic field for? (1-29):',1,29)
    else:
        day = inpInt('What day of the month do you want the magnetic field for? (1-28)',1,28)
    dayssince = 0
    if mth >1:
        for i in range(2,mth):
            dayssince = 31
            if i in day30:
                dayssince += 30
            elif i in day31:
                dayssince += 31
            elif leap == True:
                dayssince += 29
            else:
                dayssince += 28
    days = day + dayssince
    if leap == True:
        return yr + (days/366)
    else:
        return yr + (days/365)
