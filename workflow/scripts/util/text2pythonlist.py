import sys

input_file = sys.argv[1]

with open(input_file, 'r') as f:
    lines = f.readlines()

    print('lines = [')
    for line in lines:
        print(f'    "{line.rstrip()}",')
    print(']\n')

