import pandas as pd
import matplotlib.pyplot as plt
import sys




def read_data(filename):
    f = open(filename, 'r')
    data = []
    for line in f:
        data.append(line.split('\t')[:5])
    f.close()
    return data


# вызывается для кажддого файла
# возвращает словарь genes, gen - (sum of coverage, count)
def count_coverage(data):
    genes = {}
    for i in range(len(data)):
        gen = data[i][3]
        if gen not in genes.keys():
            genes[gen] = (int(data[i][4]), 1)

        else:
            old_count = genes[gen][1]
            old_reads = genes[gen][0]
            genes[gen] = (old_reads + int(data[i][4]), old_count + 1)

        if gen not in list_of_genes:
            list_of_genes.append(gen)

    return genes


# вызывается для каждого файла
# в глобальный словарь записывает среднее покрытие генома для данного файла
def count_genome_coverage(data, filename):
    sum_cov = 0
    for i in range(len(data)):
        sum_cov += int(data[i][4])

    genome_coverage[filename] = sum_cov/len(data)


# вызываетяс для кадого файла
# записывает в глоабльный словарь file и среднее покрытие этого гена в этом файле
def count_file_gen_coverage(genes, filename):
    for gen in list_of_genes:
        if gen not in file_gen_coverage.keys():
            file_gen_coverage[gen] = [(filename, genes[gen][0]/genes[gen][1])]
        else:
            file_gen_coverage[gen].append((filename, genes[gen][0] / genes[gen][1]))


# вызывается единожды
# для каждого гена подсчитывает отношение среднего покрытия гена в файле к среднему покрфтию генома файла
def count_average_gen_coverage():
    gen_coverage_by_genome = {}
    for gen in list_of_genes:
        # получили словарь
        # ген - (файл, покрытие относительно генома)
        for (file, av_cov) in file_gen_coverage[gen]:
            if gen not in gen_coverage_by_genome.keys():
                gen_coverage_by_genome[gen] = [(file, av_cov / genome_coverage[file])]
            else:
                gen_coverage_by_genome[gen].append((file, av_cov / genome_coverage[file]))

    for (file, relative_cov) in gen_coverage_by_genome['POU5F1']:
        print(file[-30:], relative_cov)


# список таргентных генов
list_of_genes = []

# среднее покрытие генома в файле
# file- av coverage genome
genome_coverage = {}

# покртие гена в файле относительно генома
# gen - [(file, sum cov/ count cov)]
file_gen_coverage = {}

# покрытие региона отнистельно генома в файле
# циферка- номер региона в бед файле-для всех одна и та же
file_region_coverage = {}


def main(filename):
    data = read_data(filename)
    genes = count_coverage(data)
    count_genome_coverage(data, filename)
    count_file_gen_coverage(genes, filename)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            main(arg)
        count_average_gen_coverage()