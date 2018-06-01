"""
Создает файл со всеми регионами для cns файла,
считывая каждый результат для каждого референса по отдельности
"""

import sys
import numpy as np

def count_common_difference(files):
    # хромосома - список циферок для каждого региона
    difference = {}
    for file in files:
        f = open(file, 'r')
        for i,line in enumerate(f):
            if i > 2 and i < 25:
                tmp = line.split('<')
                try:
                    chr = tmp[1].split(' ')[1]
                except Exception:
                    print("AAAAAA,", file)
                    print(line)
                    print(tmp)
                list01 = []
                for i in tmp[0]:
                    list01.append(float(i))
                if chr not in difference.keys():

                    difference[chr] = np.array(list01)
                else:
                    difference[chr] = difference[chr] + np.array(list01)
        f.close()


    for chr in difference.keys():
        #print(chr +' \n', np.where(difference[chr] > 5, 1, 0))
        print(chr + ' \n', difference[chr])








if __name__ == "__main__":
    if len(sys.argv) > 1:
        # первый файл - образец
        # второй - референс
        files = sys.argv[1:]
        count_common_difference(files)



