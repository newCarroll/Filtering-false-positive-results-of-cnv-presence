# смотрим cnvkit результаты
# смотрим на учаскти, которые были детектированы
# смотрим на отношение среднее покрыте региона (размера REGION_SIZE) относиельно  медианы генома в нашем и в другом образце
# если отличается - пишем в remaining cnv регион, его среднее покрытие


# вообще пока  в один регион кладу регоинычики, расстояние между двумя последовательными меньше 300




import matplotlib.pyplot as plt
import numpy as np
import sys

# четверки хромосома, начало, конец региона, ген
region_list = []


REGION_SIZE = 3000

# для каждого файла - список средних покртыий регинов
# регионы для всех файлов одинаковые
file_list_coverages = {}


# для каждого файла - список списков регионов для сооответсвующих регионов
# file_list_coverages
# регионы для всех файлов одинаковые
file_regions = {}








#TODO сделать буд файлы по нормальному бед файлу CP

#для каждой четверки - список пар - файл, покрытие
region_file_coverage_dict = {}


# среднее покрытие генома для файла
# словарь - файл: список покрытиий
file_genome_coverage = {}


# четверка - хромосома, начало и конец региона
# по идее должны идти по порядку
def create_list_regions(bedfile):
    with open(bedfile, 'r') as f:
        for line in f:
            chr, beg, end, gen = line.split('\t')[:4]
            region_list.append((chr, int(beg), int(end), gen))


# для каждого региона, для каждого файла сделаем список покрытий этого региона
def create_region_file_dict(bedfile):
    if bedfile not in file_genome_coverage.keys():
        file_genome_coverage[bedfile] = []

    count_lines = 0

    with open(bedfile, 'r') as f:
        for line in f:
            tmp = line.split('\t')
            chr, beg, end, gen = tmp[:4]
            this_file_region_coverage = int(tmp[4])
            if (chr, int(beg), int(end), gen) not in region_file_coverage_dict.keys():
                region_file_coverage_dict[(chr, int(beg), int(end), gen)] = [(bedfile, this_file_region_coverage)]
            else:
                region_file_coverage_dict[(chr, int(beg), int(end), gen)].append((bedfile, this_file_region_coverage))
            file_genome_coverage[bedfile].append(this_file_region_coverage)
            count_lines += 1

    # а если взять медиану?
    file_genome_coverage[bedfile] = np.median(np.array(file_genome_coverage[bedfile]))
    #file_genome_coverage[bedfile] = np.sum(np.array(file_genome_coverage[bedfile])) / len(file_genome_coverage[bedfile])



# нарисуем графики для каждого региона, но при этом отметим черные
# вертикальные полоски только между длинами размера REGION_SIZE,
# красные, если промежуток  больше REGION_SIZE между регионами
def print_regions_plots():
    for file in file_genome_coverage.keys():
        prev_chr = region_list[0][0]
        prev_beg, prev_end = region_list[0][1], region_list[0][2]
        # список покртыий региона, нужен,чтобы потом найти среднее  внутри него
        region_coverages = []
        # текушие регионы. размер суммы  длины регионов не превышает 300-400 пар нулеотидов
        current_regions = []

        for i,region in enumerate(region_list):
            change_chr = False
            chr, beg, end, gen = region
            coverages_list = region_file_coverage_dict[region]
            for list_file, coverage in coverages_list:
                if list_file == file:
                    current_coverage = coverage
                    break



            if chr != prev_chr:

                prev_chr = chr
                region_coverage_sum = np.sum(np.array(region_coverages))
                region_attitude = region_coverage_sum / len(region_coverages) / file_genome_coverage[file]
                if file not in file_list_coverages.keys():
                    # отношение среднего покрытие региона к среднему покрытию генома
                    file_list_coverages[file] = [region_attitude]
                    file_regions[file] = [current_regions]
                else:
                    file_list_coverages[file].append(region_attitude)
                    file_regions[file].append(current_regions)


                region_coverages = [current_coverage]
                current_regions = [region]
                prev_end = end
                prev_beg = beg
                change_chr = True


            if not change_chr:
                if beg - prev_end > REGION_SIZE:
                    region_attitude = np.sum(np.array(region_coverages)) / len(region_coverages) / file_genome_coverage[file]
                    if file not in file_list_coverages.keys():
                        # отношение среднего покрытие региона к среднему покрытию генома
                        file_list_coverages[file] = [region_attitude]
                        file_regions[file] = [current_regions]
                    else:
                        file_list_coverages[file].append(region_attitude)
                        file_regions[file].append(current_regions)

                    # если это закомментить - будет хороший результат.
                    # значит ли это что новый регион внес большой вклад ? разделить по кускам? найходить такие регионы?
                    region_coverages = [current_coverage]
                    current_regions = [region]

                else:
                    region_coverages.append(current_coverage)
                    current_regions.append(region)

                prev_end = end
                prev_beg = beg
                # если забть про это, то есть количетсво регионов только уввеличиваются -получается очень хороший результат



            #elif end > black_positions[-1] + REGION_SIZE:
            #    black_positions.append(i)




def find_gen_from_regions(list_of_strange_regions_indexes, any_file):
    # из file_regions. В нем списоег больших регионов состоящи й из спсика маленьких регионов

    for i in list_of_strange_regions_indexes:
        for (chr, beg, end, gen) in file_regions[any_file][i]:
            print(gen, end=' ')
        print()



# просто разница между геном и соседос справа сильно отличается в файлах
def find_region_neighbors(coverages1, coverages2, any_file):
    # из file_regions. В нем списоег больших регионов состоящи й из спсика маленьких регионов

    countbigger = 0
    for i in range(1, len(coverages1)-1):
        if abs((coverages1[i] / coverages1[i+1])) > 2 * abs((coverages2[i]/coverages2[i+1])):#   and abs((coverages1[i] / coverages1[i-1])) > 1.5 * abs((coverages2[i]/ coverages2[i-1])):
            #print(file_regions[any_file][i][0]) # посмотрим в разных ли они хромосомах
            countbigger += 1

    print('countbigger = ', countbigger)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        # первый файл - образец
        # второй - референс
        files = sys.argv[1:]
        for i, arg in enumerate(sys.argv[1:]):
            if i == 0:
                create_list_regions(arg)
            create_region_file_dict(arg)


        print_regions_plots()

        # получили список отношение средних покртиы й регионов относительно среднего покртыи ягенома для двух файлов
        # теперь хочу сравнить один регион с другим и посмотреть на сколько отличается

        # по регионам идущим подряд ( разделяются большим промежутком или другой хромосомой)
        rates_regions_file0 = file_list_coverages[files[0]]
        rates_regions_file1 = file_list_coverages[files[1]]

        list_strange_regions_indexes = []
        difference = np.array(rates_regions_file0) / np.array(rates_regions_file1)
        for i in range(len(difference)):
            if difference[i] >= 1.5 or difference[i] <= 0.66:
                list_strange_regions_indexes.append(i)
                #print(i, end=' ')
                #print(file_regions[files[0]][i])
        print()





        bigger = np.where(difference >= 1.5, 1, 0)
        print("biger = ", np.sum(bigger))

        smaller = np.where(difference <= 0.66, 1., 0.)
        print("smaller = ", np.sum(smaller))
        print("count of regions =",  len(difference))


        #find_gen_from_regions(list_strange_regions_indexes, files[0])
        find_region_neighbors(rates_regions_file0, rates_regions_file1, files[0])

# c ошибкой было лучше
# получились не очень хорошие показатели для двух файлов



