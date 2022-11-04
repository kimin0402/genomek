import os
import math
import subprocess
import time
from subprocess import call

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)


def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)

## Need to be checked , n==0 case
    if n == 0 :
        return 0

    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)


def getPredefinedModel(depth, Family_flag):
     if Family_flag:
         if depth > 10:
             return 0.874611,0.022596,0.644481,0.020908
         elif depth > 5:
             return 0.785312,0.021318,0.596133,0.022502
         elif depth > 2:
             return 0.650299,0.019252,0.5346,0.020694
         elif depth > 1:
             return 0.578582,0.018379,0.495017,0.021652
         elif depth > 0.5:
             return 0.524757,0.023218,0.465653,0.027378
         else:
    #         print "Warning: Sample region depth is too low < 1"
             return 0.524757,0.023218, 0.465653, 0.027378
     else:
         if depth > 10:
             return 0.874546, 0.022211, 0.310549, 0.060058
         elif depth > 5:
             return 0.785249,0.021017, 0.279778, 0.054104
         elif depth > 2:
             return 0.650573, 0.018699,0.238972, 0.047196
         elif depth > 1:
             return 0.578386,0.018526, 0.222322, 0.041186
         elif depth > 0.5:
             return 0.529327,0.025785, 0.217839, 0.040334
         else:
    #         print "Warning: Sample region depth is too low < 1"
             return 0.529327,0.025785, 0.217839, 0.040334


def classifyNV(vec2Classify, p0Vec, p0S, p1Vec, p1S):
    if abs(p0Vec - vec2Classify) - p0S > abs(p1Vec - vec2Classify) - p1S:
        return abs((abs(p0Vec - vec2Classify) - p0S )/ (abs(p1Vec - vec2Classify) -  p1S )), 1
    else: 
        return abs((abs(p0Vec - vec2Classify) - p0S) / (abs(p1Vec - vec2Classify)  -  p1S)), 0  


def run_ncm(vcf_list, ref_ver='37', outdir='test/'):
    '''
    vcf_list: list of vcf paths
    '''
    glob_scores = dict()    #Whole score
    feature_list = dict()   #Each Feature List
    label = []              #Samples
    features = []           #dbSNP features
    mean_depth = dict()
    real_depth = dict()
    real_count = dict()
    sum_file = dict()
    out_tag = ""
    pdf_tag = ""
    Family_flag = False
    Nonzero_flag = False

    if ref_ver == '37':
        bedFile = "/home/users/kimin/tools/NGSCheckMate/SNP/SNP_GRCh37_hg19_woChr.bed"
        # bedFile = "/home/users/kimin/tools/NGSCheckMate/SNP/SNP_GRCh37_hg19_wChr.bed"
    else:
        bedFile = "/home/users/kimin/tools/NGSCheckMate/SNP/SNP_GRCh38_hg38_wChr.bed"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    key_order = open(bedFile,'r')
    for line in key_order.readlines():
        if line.startswith("#"):
            continue
        temp = line.strip().split('\t')
        if temp[0].find("chr")!= -1:
            features.append(str(temp[0][3:])+"_"+ str(temp[2]))
        else:   
            features.append(str(temp[0])+"_"+ str(temp[2]))

    ## adopted from function 'createDataSetFromList'
    for link in vcf_list:
        f = open(link, "r")
        dbsnpf= open(bedFile,"r")
        file = link[link.rindex("/")+1:]
        depth = dict()
        depth[file] = 0
        real_count[file] = 0
        count = 0

        sum_dict=dict()
        sum_dict[file] = 0
        scores = dict()     # Scores of B-allel Frequencies
        #DBSNP ID collecting system
        for i in dbsnpf.readlines():
            temp = i.strip().split('\t')
            if temp[0].find("chr")!= -1:
                ID = str(temp[0][3:]) + "_" + str(temp[2])
            else:   
                ID = str(temp[0]) + "_" + str(temp[2])
            scores[ID] = 0
            count = count + 1

        ## 0618_samtools and haplotyper
        vcf_flag = 0

        score_set = dict()
        #VCF file PROCESSING  and Generation of features
        total = 0
        GVCF_samples = dict()
        for i in f.readlines():        
            if i.startswith("#"):
                if i.find("DP4") != -1:
                    vcf_flag = 1
                if i.find("#CHROM") != -1:
                    temp = i.strip().split('\t')
                    total=len(temp) - 9
                    if total != 1:
                        for sample_idx in range(0,total):
                            file = temp[sample_idx + 9]
                            GVCF_samples[temp[sample_idx + 9]] = []
                            score_set[temp[sample_idx + 9]] = dict()
                            depth[temp[sample_idx + 9]] = 0
                            real_count[temp[sample_idx + 9]] = 0
                            sum_dict[temp[sample_idx + 9]] =0
                            feature_list[temp[sample_idx + 9]] = []
                    if total == 1:
                        feature_list[file] = []
                continue

            temp = i.strip().split('\t')
          ## ID in BED file only 
            if temp[0].find("chr")!= -1:
                ID = str(temp[0][3:]) + "_" + str(temp[1])
            else:   
                ID = str(temp[0]) + "_" + str(temp[1])

            if ID not in scores:
                continue

            if vcf_flag == 1:
                values = temp[7].split(';')

                if values[0].startswith("INDEL"):
                    continue

                for j in values:
                    if j.startswith("DP4"):
                        readcounts = j.split(',')
                        readcounts[0] = readcounts[0][4:]
                        total_reads =(float(readcounts[0]) + float(readcounts[1]) + float(readcounts[2]) + float(readcounts[3])) 
                        score = 0
                        if total_reads > 0:
                            score = (float(readcounts[2]) + float(readcounts[3])) / total_reads
                            real_count[file] = real_count[file] + 1  

                        depth[file] = depth[file] + total_reads

                        if ID in scores:
                            feature_list[file].append(ID)
                            scores[ID]= score
                            sum_dict[file] = sum_dict[file] + float(readcounts[2]) + float(readcounts[3])
            elif total == 1 and vcf_flag == 0:
                format = temp[8].split(':')  ##Format
                AD_idx = -1 
                DP_idx = -1
                for idx in range(0,len(format)):
                    if format[idx] == "AD":
                        AD_idx = idx
                    elif format[idx] == "DP":
                        DP_idx = idx
                if AD_idx == -1:
                    continue
                if DP_idx == -1:
                    continue
                idx = 9
                values = temp[idx].split(":")
                readcounts = values[AD_idx].split(',')

                if float(readcounts[0]) + float(readcounts[1]) < 0.5:
                    score =0
                else:
                    score = float(readcounts[1])/ (float(readcounts[0]) + float(readcounts[1]))
                depth[file] = depth[file] + float(values[DP_idx])
                if float(values[DP_idx]) > 0:
                    real_count[file] = real_count[file] + 1
                if ID in scores:
                    feature_list[file].append(ID)
                    scores[ID]= score  ##from here!
                    sum_dict[file] = sum_dict[file] + float(readcounts[1])
            else:  ###### Haplotyper or other VCF
                format = temp[8].split(':')  ##Format
                AD_idx = -1
                DP_idx = -1
                for idx in range(0,len(format)):
                    if format[idx] == "AD":
                        AD_idx = idx
                    elif format[idx] == "DP":
                        DP_idx = idx
                if AD_idx == -1:
                    continue
                if DP_idx == -1:
                    continue
                idx = 9
                for file in GVCF_samples:
                    values = temp[idx].split(":")
                    if len(values) < len(format):
                        score = 0
                        idx = idx + 1
                        continue

                    readcounts = values[AD_idx].split(',')

                    if float(readcounts[0]) + float(readcounts[1]) < 0.5:
                        score =0
                    else:
                        score = float(readcounts[1])/ (float(readcounts[0]) + float(readcounts[1]))
                    depth[file] = depth[file] + float(values[DP_idx])
                    if float(values[DP_idx]) > 0:
                        real_count[file] = real_count[file] + 1                   

                    if ID in scores:
                        feature_list[file].append(ID)
                        score_set[file][ID]= score   ##from here!
                        sum_dict[file] = sum_dict[file] + float(readcounts[1])

                    idx = idx + 1
                    
        if total == 1:
            mean_depth[file] = depth[file] / float(count)
            real_depth[file] = depth[file] / float(real_count[file])
            sum_file[file] = sum_dict[file]

            for key in features:
                if file in glob_scores:
                    glob_scores[file].append(scores[key])
                else:
                    glob_scores[file] = [scores[key]]    
        else:
            for file in GVCF_samples:
                mean_depth[file] = depth[file] / float(count)
                real_depth[file] = depth[file] / float(real_count[file])
                sum_file[file] = sum_dict[file]

                for key in features:
                    if key not in score_set[file]:
                        score_set[file][key] = 0
                    if file in glob_scores:
                        glob_scores[file].append(score_set[file][key])
                    else:
                        glob_scores[file] = [score_set[file][key]]    
        dbsnpf.close()
        f.close()

    for key in sorted(glob_scores):
        label.append(key)

    ## adopted from function classifying
    AUCs =[]
    wholeFeatures = 50
    temp =[]
    altFreqList = []
    keyList = []

    for key in sorted(glob_scores):
        altFreqList.append(glob_scores[key])
        keyList.append(key)

    dataSetSize = len(altFreqList)
    filter_list = []

    for i in range(0, dataSetSize):
        for j in range(0, dataSetSize):
            if i!=j:
                if keyList[j] not in filter_list:
                    temp.append([keyList[i],keyList[j]])
        filter_list.append(keyList[i])

    for iterations in range(49,wholeFeatures):

        samples = []
        numFeatures = iterations
        count = 0

        for i in range(0,len(temp)):
            tempA = set(feature_list[temp[i][0].strip()])
            tempB = set(feature_list[temp[i][1].strip()])

            selected_feature = tempA.intersection(tempB)
            vecA = []
            vecB = []

            idx = 0
            for k in features:
                if k in selected_feature:
                    vecA.append(glob_scores[temp[i][0].strip()][idx])
                    vecB.append(glob_scores[temp[i][1].strip()][idx])
                idx = idx + 1

            distance = pearson_def(vecA, vecB)
            samples.append(distance)

        predStrength = []
        training_flag =0
    ####0715 Append

        output_matrix_f = open(outdir + "/output_corr_matrix.txt","w")
        output_matrix = dict()

        if out_tag!="stdout":
            out_f = open(outdir + "/" + out_tag + "_all.txt","w")
            out_matched = open(outdir + "/" + out_tag + "_matched.txt","w")

        for i in range(0, len(keyList)):
            output_matrix[keyList[i]] = dict()
            for j in range(0,len(keyList)):
                output_matrix[keyList[i]][keyList[j]] = 0

        if training_flag == 1:
            #make training set
            for i in range(0,len(samples)):
                trainMatrix= []
                trainCategory = []
                for j in range(0, len(samples)):
                    if i==j:
                        continue
                    else:
                        trainMatrix.append(samples[j])
                        trainCategory.append(classLabel[j])
                #training samples in temp
                #p0V, p1V, pAb = trainNB0(array(trainMatrix),array(trainCategory))
                p1V,p1S, p0V, p0S = trainNV(array(trainMatrix),array(trainCategory))
                result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
                if result[1] == 1:
                    print(str(temp[i][0]) + '\tsample is matched to\t',str(temp[i][1]),'\t', samples[i])
                predStrength.append(result[0])
        else :
            for i in range(0,len(samples)):
                depth = 0
                if Nonzero_flag: 
                    depth = min(real_depth[temp[i][0].strip()],real_depth[temp[i][1].strip()])
                else:
                    depth = min(mean_depth[temp[i][0].strip()],mean_depth[temp[i][1].strip()])

                p1V,p1S, p0V, p0S = getPredefinedModel(depth, Family_flag)
                result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
                if result[1] ==1:
                    output_matrix[temp[i][0].strip()][temp[i][1].strip()] = samples[i]
                    output_matrix[temp[i][1].strip()][temp[i][0].strip()] = samples[i]
                    if out_tag=="stdout":
                        print(str(temp[i][0]) + '\tmatched\t',str(temp[i][1]),'\t', round(samples[i],4),'\t',round(depth,2))
                    else :
                        out_f.write(str(temp[i][0]) + '\tmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                        out_matched.write(str(temp[i][0]) + '\tmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                else:
                    if out_tag=="stdout":
                        print(str(temp[i][0]) + '\tunmatched\t',str(temp[i][1]),'\t', round(samples[i],4),'\t',round(depth,2))
                    else :
                        out_f.write(str(temp[i][0]) + '\tunmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                predStrength.append(result[0])
            #testing sample is samples
        output_matrix_f.write("sample_ID")
        for key in output_matrix.keys():
            if key.find(".vcf") != -1:
                output_matrix_f.write("\t" + key[0:key.index('.vcf')])
            else:
                output_matrix_f.write("\t" + key)
        output_matrix_f.write("\n")

    #        for key in output_matrix.keys():
    #            for otherkey in output_matrix[key].keys():
    #                if output_matrix[key][otherkey] != 0:
    #                    output_matrix[otherkey][key] = output_matrix[key][otherkey] 

        for key in output_matrix.keys():
            if key.find(".vcf") != -1:
                output_matrix_f.write(key[0:key.index('.vcf')])
            else:
                output_matrix_f.write(key)
            for otherkey in output_matrix.keys():
                output_matrix_f.write("\t" + str(output_matrix[key][otherkey]))
            output_matrix_f.write("\n")   

        output_matrix_f.close()         
        if out_tag!="stdout":
            out_f.close()
            out_matched.close()   