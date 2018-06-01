# для каждого вида нозологии для каждого гена запоминаем вероятность
# амплификации и делеции

import math
import numpy as np



# файл - словарь для генов
file_dict = {}
def analyse_gen_probability(filename):
    f = open("oncotree/" + filename, 'r')
    cancer_kind = filename.split('.')[0].split('_')[-1]
    file_dict[cancer_kind] = {}
    for i, line in enumerate(f):
        if i > 0:
            tmp = line.split('\t')
            gen = tmp[0]
            # такая маленькая? /100
            # в письме и файле разная информация где делеция где амлификация
            del_prob = (float(tmp[1]) + float(tmp[2])) / 2 / 100
            gain_prob = (float(tmp[3]) + float(tmp[4])) / 2 / 100
            file_dict[cancer_kind][gen] = (del_prob, gain_prob)
    f.close()


# задается список файлов референсов
def count_common_difference_deletion(files):
    # хромосома,номерок bigregion  - список циферок для каждого региона
    # цифра в регионе - количество файлов, где отличается покрытие в два и более раза в нижнюю сторону
    difference = {}
    for file in files:
        f = open('accross_cnvkit/results/' + file, 'r')
        for i,line in enumerate(f):
            if i > 2 and i < 72:
                tmp = line.split(' ')
                chr = tmp[2]
                list01 = []
                for j in tmp[0]:
                    list01.append(float(j))
                if (chr,i) not in difference.keys():
                    difference[(chr,i)] = np.array(list01)
                else:
                    difference[(chr,i)] = difference[(chr,i)] + np.array(list01)
        f.close()

    return difference


def create_list_genes():
    gen_list = []
    f = open('oncotree/_max_LIVER.txt', 'r')
    for i, line in enumerate(f):
        if i > 0:
            tmp = line.split('\t')
            gen_list.append(tmp[0])
    f.close()
    return gen_list



# вероятность такого отностельнго покртия региона при сnv-delteion
# 10 при cnv - совсем маленькая
# 9 чуть больше
# лямбда - среднее отношение покрытия регионов отноистельно хромосом

# data nocnv
# при 10 большая
# при 9 чуть меньше
# при 0.5 совсем маленькая




# правда ли, что мы смотрим увеличение отноистельно начала хромосомы
def model_rate_coverage_amplification(region_coverage, gen, cancer_kind):
    # задать lambda_
    # относительное покрытие региона в хромосоме
    if region_coverage != 0:
        k = round(region_coverage)
        # стоит ли вообще округлять
        region_probability_data_nocnv = 1**k / math.factorial(k) * math.exp(-1)
        # region_probability_data_cnv = 10 **k / math.factorial(k) * math.exp(-10)
        region_probability_data_cnv = 1 - region_probability_data_nocnv
        p_cnv = file_dict[cancer_kind][gen][0]
        p_nocnv = 1 - p_cnv
        rate_deletion = region_probability_data_cnv * p_cnv / (region_probability_data_nocnv * p_nocnv)
        return rate_deletion
    else:
        return -1


# правда ли, что мы смотрим увеличение отноистельно начала хромосомы
def model_rate_coverage_deletion(region_coverage, gen, cancer_kind):
    # задать lambda_
    lambda_ = 1
    # относительное покрытие региона в хромосоме
    if region_coverage != 0:
        k = round(1 / region_coverage)
        # стоит ли вообще округлять
        if k > 20:
            k = 20
        region_probability_data_nocnv = lambda_ ** k / math.factorial(k) * math.exp(-lambda_)
        region_probability_data_cnv = 1 - region_probability_data_nocnv
        p_cnv = file_dict[cancer_kind][gen][0]
        p_nocnv = 1 - p_cnv
        rate_deletion = region_probability_data_cnv * p_cnv / region_probability_data_nocnv * p_nocnv
        return rate_deletion
    else:
        return -1



regions_cnvkit = []
bigregions = []
bigregion_coverages = []
# регионы cnvkit - Характеризуются тройкой (хромосома, начало, конец)
def create_cnvkit_regions_list(cnvkit_file):
    f = open(cnvkit_file, 'r')
    for i, line in enumerate(f):
        if i > 0:
            tmp = line.split('\t')
            chr, beg, end = tmp[0], float(tmp[1]), float(tmp[2])
            regions_cnvkit.append((chr, beg, end))
            bigregion_coverages.append([])
            bigregions.append([])
    f.close()


region_list = []
# четверка - хромосома, начало и конец региона, ген
# по идее должны идти по порядку
def create_list_regions(bedfile):
    with open(bedfile, 'r') as f:
        for line in f:
            chr, beg, end, gen = line.split('\t')[:4]
            region_list.append((chr, int(beg), int(end), gen))


region_sample_coverage_dict = {}
# для каждой четверки==региона в базовом файле покрытие
def create_region_sample_coverage_dict(sample_bedfile):
    with open(sample_bedfile, 'r') as f:
        for line in f:
            tmp = line.split('\t')
            chr, beg, end, gen = tmp[:4]
            this_file_region_coverage = int(tmp[4])
            region_sample_coverage_dict[(chr, int(beg), int(end), gen)] = this_file_region_coverage


REGION_SIZE = 300
# посчитаем среднее для регионов файла, размер которых примерно равен REGION_SIZE, и расстояние меньше REGION_SIZE
def count_region_size_region_coverages():

    prev_chr = region_list[0][0]
    prev_beg, prev_end = region_list[0][1], region_list[0][2]
    current_region_coverages = []
    current_regions = []
    length_sum_region = 0

    region_size_regions = []
    region_size_regions_coverages = []

    for i, region in enumerate(region_list):
        change_chr = False
        chr, beg, end, gen = region
        current_coverage = region_sample_coverage_dict[region]

        if chr != prev_chr:

            prev_chr = chr
            region_coverage_sum = np.sum(np.array(current_region_coverages))
            # среднее покрытие маленького региона
            region_attitude = region_coverage_sum / len(current_region_coverages)
            region_size_regions_coverages.append(region_attitude)
            region_size_regions.append(current_regions)

            current_region_coverages = [current_coverage]
            current_regions = [region]
            prev_end = end
            prev_beg = beg
            change_chr = True
            length_sum_region = end - beg

        if not change_chr:
            if beg - prev_end > REGION_SIZE or length_sum_region > REGION_SIZE:
                region_attitude = np.sum(np.array(current_region_coverages)) / len(current_region_coverages)
                region_size_regions_coverages.append(region_attitude)
                region_size_regions.append(current_regions)
                current_region_coverages = [current_coverage]
                current_regions = [region]
                length_sum_region = end - beg
            else:
                current_region_coverages.append(current_coverage)
                current_regions.append(region)
                length_sum_region += end - beg

            prev_end = end
            prev_beg = beg

    region_attitude = np.sum(np.array(current_region_coverages)) / len(current_region_coverages)
    region_size_regions_coverages.append(region_attitude)
    region_size_regions.append(current_regions)

    # print(region_size_regions[6:10])
    # print(region_size_regions_coverages[6:10])

    for i, cnvkit_region in enumerate(regions_cnvkit):
        cnk_chr, cnk_beg, cnk_end = cnvkit_region
        for j, some_regions in enumerate(region_size_regions):
            reg_chr, reg_b, reg_end, reg_gen = some_regions[0]
            if reg_chr == cnk_chr:
                if reg_b >= cnk_beg and reg_end <= cnk_end:
                    bigregions[i].append(some_regions)
                    bigregion_coverages[i].append(region_size_regions_coverages[j])

    # print(bigregions[18][8:])
    # print(bigregion_coverages[18][8:])


def create_regions_list(sample_bed_file, cns_file):
    create_cnvkit_regions_list(cns_file)
    create_list_regions(sample_bed_file)
    create_region_sample_coverage_dict(sample_bed_file)
    count_region_size_region_coverages()



def find_part_chrom_coverage(part):
    chrom_part_dict = {}
    prev_chr = region_list[0][0]
    current_region_coverages = []

    for i, region in enumerate(region_list):
        change_chr = False
        chr, beg, end, gen = region
        current_coverage = region_sample_coverage_dict[region]

        if chr != prev_chr:
            l = len(current_region_coverages)
            chrom_part_coverage = np.mean(np.array(current_region_coverages[:int(l * part)]))
            chrom_part_dict[prev_chr] = chrom_part_coverage

            current_region_coverages = [current_coverage]
            change_chr = True
            prev_chr = chr

        if not change_chr:
            current_region_coverages.append(current_coverage)

    l = len(current_region_coverages)
    chrom_part_coverage = np.mean(np.array(current_region_coverages[:int(l * part)]))
    chrom_part_dict[prev_chr] = chrom_part_coverage
    return chrom_part_dict


def main():
    file_list = ["_max_ADRENAL_GLAND.txt","_max_AMPULLA_OF_VATER.txt","_max_BILIARY_TRACT.txt","_max_BLADDER.txt",
                 "_max_BLOOD.txt","_max_BONE.txt","_max_BOWEL.txt","_max_BRAIN.txt","_max_BREAST.txt","_max_CERVIX.txt",
                 "_max_EYE.txt","_max_HEAD_NECK.txt","_max_KIDNEY.txt","_max_LIVER.txt","_max_LUNG.txt",
                 "_max_LYMPH.txt","_max_OTHER.txt","_max_OVARY.txt","_max_PANCREAS.txt","_max_PENIS.txt",
                 "_max_PERITONEUM.txt","_max_PLEURA.txt","_max_PNS.txt","_max_PROSTATE.txt","_max_SKIN.txt",
                 "_max_SOFT_TISSUE.txt","_max_STOMACH.txt","_max_TESTIS.txt","_max_THYMUS.txt","_max_THYROID.txt",
                 "_max_UTERUS.txt","_max_VULVA.txt"]

    for filename in file_list:
        analyse_gen_probability(filename)

#-----------------------------------------------------------------------------------------------------------------------


    # ПОЛУЧИМ РЕГИОНЫ СООТВЕТВУЮЩИЕ РЕГИОНАМ cnvkit
    sample_bed_file = 'bed_files/sample_7/coverage_mil-2017-ccp-mil_k3_7_pcbl-t.sorted.bam.bed'
    cns_file = 'mil-2017-ccp-mil_k3_7_pcbl-t.sorted.call.cns'
    # список регионов соотвествующих cnvkit, каждый состоит из регионов длины REGION_SIZE
    # соответвующие списки средних покрытий
    # пока ищем для всех регионов
    # в том числе не увеличенных\уменьшенных
    create_regions_list(sample_bed_file, cns_file)



    # ПОЛУЧИМ СРЕДНЕЕ ПОКРЫТИЕ ДЛЯ КАЖДОЙ ХРОМОСОМЫ ДЛЯ ЕЕ ПЕРВОЙ Part части
    part = 1
    chrom_part_dict = find_part_chrom_coverage(part)


#-----------------------------------------------------------------------------------------------------------------------
    # НАШЛИ ВЕРОЯТНОСТИ ДЕЛЕЦИИ ДЛЯ КАЖДОГО РЕГОНЧИКА В cnvkit РЕГИОНЕ
    # bigregions выглядит так примерно -  [[300 [[][][]], 300 [[][]]] - cnvregion]
    # для каждого региона размера 300 знаем его среднее покрытие
    cancer_kind = 'TISSUE'
    rate_probabilites = [[] for i in range(len(bigregions))]
    for i, bigregion in enumerate(bigregions[:-1]):
        for j, region in enumerate(bigregion):  # <--- проходим по регионам размера 300
            chr = region[0][0]
            gen = region[0][3]
            region_coverage = bigregion_coverages[i][j] / chrom_part_dict[chr]  # <--- посмотрели относительно начала хромосомы
            if region_coverage > 1:
                # if gen in file_dict[cancer_kind].keys():
                #     region_prob_rate_amplification = model_rate_coverage_amplification(region_coverage, gen, cancer_kind)
                #     #region_prob_rate_amplification = model_rate_coverage_deletion(region)
                # else:
                region_prob_rate_amlification = '_'
            else:
                if gen in file_dict[cancer_kind].keys():
                    region_prob_rate_deletion = round(model_rate_coverage_deletion(region_coverage, gen, cancer_kind),4)
                else:
                    region_prob_rate_deletion = '_'

            #region_prob_rate_amplification = model_rate_coverage_amlification(region)
            rate_probabilites[i].append(region_prob_rate_deletion)
            #print(round(region_prob_rate_deletion, 2), end=' ')
            print(region_prob_rate_deletion, end=' ')

        print('<---', regions_cnvkit[i][0])






main()