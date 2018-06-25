#!/usr/bin/env python3

def main():
    test1 = {123:['Tr1','Tr2','Tr3'], 45:['Tr4','Tr5']}
    result = {}
    for k,v in test1.items():
        for value in v:
            result.update({value:k})
    print(result)

if __name__ == '__main__':
    main()