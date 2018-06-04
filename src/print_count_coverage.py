import pandas as pd
import matplotlib.pyplot as plt
import sys


"""
скрипт ищет гены, который сильно отличются по покрытию от среднего в геноме,
а именно выбросы
создает csv файл с такими генами - strange genes
"""


genes_coverage = {}
# cписки покрытий для генов
# обнуляются для каждого файла


def read_data(filename):
    f = open(filename, 'r')
    data = []
    for line in f:
        data.append(line.split('\t')[:5])
    f.close()
    return data


def count_coverage(data):
    genes = {}
    chrom = {}
    for i in range(len(data)):
        gen = data[i][3]
        chr = data[i][0]
        if gen not in genes.keys():
            genes[gen] = (int(data[i][4]), 1)
            if chr not in chrom:
                chrom[chr] = [gen]
            else:
                chrom[chr].append(gen)
            genes_coverage[gen] = [int(data[i][4])]

        else:
            old_count = genes[gen][1]
            old_reads = genes[gen][0]
            genes[gen] = (old_reads + int(data[i][4]), old_count + 1)
            genes_coverage[gen].append(int(data[i][4]))

    return genes, chrom


# принмаиет список покртыий по порядку для гена, ген, хромосому
def gen_plot(list_gen_coverage, gen, chr, filename):
    x = [i for i in range(len(list_gen_coverage))]
    l = len(list_gen_coverage)
    sorted_gen_coverage = sorted(list_gen_coverage)
    label = gen+' '+chr
    plt.scatter(x, list_gen_coverage, c='gray', s=2, label=label)
    plt.axhline(sorted_gen_coverage[int(len(list_gen_coverage)/2)], c='green')
    plt.legend()
    x_positions = [i*10 for i in range(int(len(x)/10))]
    for x_position in x_positions:
        plt.axvline(x_position)

    index_q3 = int(l / 4 * 3)
    index_q1 = int(l / 4)
    q3 = sorted_gen_coverage[index_q3]
    q1 = sorted_gen_coverage[index_q1]

    diap = q3 - q1
    low_border = diap * 1.5 + q3
    high_border = diap * 3 + q3

    down_light_border = q1 - diap * 1.5
    down_strong_border = q1 - 3 * diap
    plt.axhline(high_border, c='red')
    # нижняя граница
    plt.axhline(down_strong_border, c='red')
    # считается, что cnv от 1000 нуклеотидов
    # один регион примерно 100 нуклеотидов
    # нарисуем вертикальный линии на каждых 1000 - то есть примерно каждый десятый
    picture_name = filename + '.plots/' + chr + '_' + gen + '.png'
    plt.savefig(picture_name)
    plt.clf()


def main(filename):
    print(filename)
    cov8 = read_data(filename)
    genes , chrom = count_coverage(cov8)


    for chr in chrom.keys():
        chr_genes = chrom[chr]
        cov = []
        for gen in chr_genes:
            cov.append((genes[gen][0]/genes[gen][1], gen))

        cov = sorted(cov, key=lambda key: key[0])
        #cov = sorted(cov)
        l = len(cov)


        x_positions = []
        # список позиций где гены заканыиватся
        # хочу отделить их друг от друга


        # список покрытий для генов в хромосоме
        #  по порядку следования генов в бед файле,
        list_of_coverages = []
        for gen in genes:
            if gen in chrom[chr]:
                list_of_coverages.extend(genes_coverage[gen])
                x_positions.append(len(list_of_coverages))
                #gen_plot(genes_coverage[gen], gen, chr, filename)

        x = [i for i in range(len(list_of_coverages))]

        # сам график покрытий
        plt.scatter(x, list_of_coverages, s=2, c='silver', label=chr)
        plt.legend()
        # медиана
        plt.axhline(cov[int(len(cov)/2)][0], c='green')
        # разделеяем гены
        for x_position in x_positions:
            plt.axvline(x_position, c='gray')

        index_q3 = int(l/4*3)
        index_q1 = int(l/4)
        q3 = cov[index_q3][0]
        q1 = cov[index_q1][0]

        diap = q3-q1
        low_border = diap * 1.5 + q3
        high_border = diap * 3 + q3

        down_light_border = q1 - diap *1.5
        down_strong_border = q1 - 3*diap

        # то есть если есть гены покртие которых сильно отличается
        #  - мыузнаем об этом - для каждого гена

        chrom_strange_genes = []
        print('\n' + chr, end=' ')
        for (c, gen) in cov:
            if c > high_border:
                chrom_strange_genes.append('so' + gen)
            elif c > low_border:
                chrom_strange_genes.append(gen)

            if c < down_strong_border:
                chrom_strange_genes.append('so' + gen)
            elif c < down_strong_border:
                chrom_strange_genes.append('so' + gen)

        if filename not in chrom_strange_genes_dict.keys():
            chrom_strange_genes_dict[filename] = [(chr, chrom_strange_genes)]
        else:
            chrom_strange_genes_dict[filename].append((chr, chrom_strange_genes))

        """
        # верхняя граниица
        plt.axhline(high_border, c='red')
        # нижняя граница
        plt.axhline(down_strong_border, c='red')
        picture_name = filename + '.plots/' + chr +'.png'
        plt.savefig(picture_name)
        plt.clf()

        #plt.show()
        """

    print()


filename_coverage = {}

# словарь, где каждого файла для каждой хромосоме будет выведен список
# выставляющихся генов, с пометкой so если сильно выставляются
chrom_strange_genes_dict = {}

chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
               "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX" ]



if __name__ == '__main__':
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            genes_coverage = {}
            main(arg)
            
        # добавим недостающих хромосом для файла

        for filename in sys.argv[1:]:
            file_existed_chr = []
            for (chr, genes) in chrom_strange_genes_dict[filename]:
                file_existed_chr.append(chr)
            for chr in chromosomes:
                if chr not in file_existed_chr:
                    chrom_strange_genes_dict[filename].append((chr, []))

        

        w = open('strange_genes_by_chrom_and_files.csv', 'w+')
        w.write('sample\t')
        for chr in chromosomes:
            w.write(chr +'\t')
        w.write('\n')
        for filename in sys.argv[1:]:
            w.write(filename+'\t')

            sorted_array = []
            for chr in chromosomes:
                for (file_chr, gen_list) in chrom_strange_genes_dict[filename]:
                    if file_chr == chr:
                        sorted_array.append((file_chr, gen_list))

            for (chr, genes) in sorted_array:
                for gen in genes:
                    w.write(gen + ',')
                w.write('\t')
            w.write('\n')

        w.close()