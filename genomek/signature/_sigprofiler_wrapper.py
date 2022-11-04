import os
import sys
import tempfile
import shutil
import subprocess
import collections
import itertools
import pandas as pd
import numpy as np
import pickle as pkl
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matgen   
from sigProfilerPlotting import sigProfilerPlotting as sigPlt_original
from sigProfilerPlotting import sample_portrait as sigPlt_portrait
from SigProfilerExtractor import subroutines as sub
from SigProfilerExtractor import single_sample as ss
from SigProfilerExtractor import PlotDecomposition as sp
from SigProfilerExtractor import SigProfilerPlottingMatrix as sigPlt
import SigProfilerExtractor as cosmic
from sigproSS import spss 
from PyPDF2 import PdfFileMerger
from PIL import Image, ImageDraw, ImageFont
from pdf2image import convert_from_path


paths = cosmic.__path__[0]






#########################
## SigProfiler Wrapper ##
#########################


def data_generation(vcf_path):
    '''
    Wrapper function of python sigprofiler package for interactive use.
    This function takes a vcf path as an input and returns a data dictionary for further use.
    '''
    
    ## Create a temporary folder and move the input vcf file to the folder
    ## Decompress using bgzip 
    ## Run matrix generation and add some additional keys
    with tempfile.TemporaryDirectory() as tempdir:
        print("Temporary directory created: ", tempdir)
        shutil.copy(vcf_path, tempdir)
        filename = tempdir+"/"+os.path.basename(vcf_path)
        cmd = f'bgzip -d {filename}'    
        subprocess.run(cmd, shell=True)
        data = matgen.SigProfilerMatrixGeneratorFunc('mat', 'GRCh37', tempdir)
        filename = tempdir+"/output/DBS/"+"mat.DBS186.all"
        try:
            data['DBS186'] = pd.read_csv(filename, sep='\t', index_col="MutationType")
        except:
            data['DBS186'] = None
        filename = tempdir+"/output/ID/"+"mat.ID28.all"
        try:
            data['ID28'] = pd.read_csv(filename, sep='\t', index_col="MutationType")
        except:
            data['ID28'] = None
        filename = tempdir+"/output/ID/"+"mat.ID415.all"
        try:
            data['ID415'] = pd.read_csv(filename, sep='\t', index_col="MutationType")
        except:
            data['ID415'] = None

    print("Dictionary is returned with following keys: ", data.keys())

    return data


def plot_sig(sample_DF, percentage=False):
    '''
    Given a dataframe, a signature plot of columns is drawn.
    '''
    plot_is_dictionary = False
    plot_is_pdf = False
    mtype = str(sample_DF.shape[0])
    sample_DF_reset = sample_DF.reset_index()

    with tempfile.TemporaryDirectory() as tempdir:
        if mtype in ['96', '192', '288', '384', '1536']:
            plot_dict = sigPlt.plotSBS(sample_DF_reset, tempdir, 'tmp', mtype, percentage)
            plot_is_dictionary = True
        elif mtype in ['78']:
            plot_dict = sigPlt.plotDBS(sample_DF_reset, tempdir, 'tmp', mtype, percentage)
            plot_is_dictionary = True
        elif mtype in ['83']:
            plot_dict = sigPlt.plotID(sample_DF_reset, tempdir, 'tmp', mtype, percentage)
            plot_is_dictionary = True

        elif mtype in ['6']:
            df2csv(df=sample_DF, fname=f'{tempdir}/6.tsv')
            sigPlt_original.plotSBS(f'{tempdir}/6.tsv', f'{tempdir}/', 'tmp', mtype, percentage)
            plot_pdf = convert_from_path(f'{tempdir}/SBS_6_plots_tmp.pdf', dpi=200)
            plot_is_pdf = True
        elif mtype in ['12', '24']:
            df2csv(df=sample_DF, fname=f'{tempdir}/24.tsv')
            sigPlt_original.plotSBS(f'{tempdir}/24.tsv', f'{tempdir}/', 'tmp', mtype, percentage)
            plot_pdf = convert_from_path(f'{tempdir}/SBS_24_plots_tmp.pdf', dpi=200)
            plot_is_pdf = True
        elif mtype in ['4608']:
            df2csv(df=sample_DF, fname=f'{tempdir}/4608.tsv')
            sigPlt_original.plotSBS(f'{tempdir}/4608.tsv', f'{tempdir}/', 'tmp', mtype, percentage)
            plot_pdf = convert_from_path(f'{tempdir}/SBS_4608_plots_tmp.pdf', dpi=200)
            plot_is_pdf = True

        elif mtype in ['94']:
            df2csv(df=sample_DF, fname=f'{tempdir}/83.tsv')
            sigPlt_original.plotID(f'{tempdir}/83.tsv', f'{tempdir}/', 'tmp', mtype, percentage)
            plot_pdf = convert_from_path(f'{tempdir}/ID_83_plots_tmp.pdf', dpi=200)
            plot_is_pdf = True
        elif mtype in ['28']:
            df2csv(df=sample_DF, fname=f'{tempdir}/28.tsv')
            sigPlt_original.plotID(f'{tempdir}/28.tsv', f'{tempdir}/', 'tmp', mtype, percentage)
            plot_pdf = convert_from_path(f'{tempdir}/ID_simple_plots_tmp.pdf', dpi=200)
            plot_is_pdf = True
        elif mtype in ['415']:
            df2csv(df=sample_DF, fname=f'{tempdir}/415.tsv')
            sigPlt_original.plotID(f'{tempdir}/415.tsv', f'{tempdir}/', 'tmp', mtype, percentage)
            plot_pdf = convert_from_path(f'{tempdir}/ID_TSB_plots_tmp.pdf', dpi=200)
            plot_is_pdf = True

        elif mtype in ['186', '312']:
            df2csv(df=sample_DF, fname=f'{tempdir}/186.tsv')
            sigPlt_original.plotDBS(f'{tempdir}/186.tsv', f'{tempdir}/', 'tmp', mtype, percentage)
            plot_pdf = convert_from_path(f'{tempdir}/DBS_186_plots_tmp.pdf', dpi=200)
            plot_is_pdf = True

    if plot_is_dictionary:
        plot_pngs = [Image.open(v) for k, v in plot_dict.items()]
        plot_pngs_combined = Image.fromarray(np.vstack((np.asarray(i) for i in plot_pngs)))
        print("Returned object is a PIL image.")
        return plot_pngs_combined

    if plot_is_pdf:
        print("Returned object is a list. Index it to show plots.")
        return plot_pdf

    else:
        print(f"Mutation type: {mtype}")
        print("This mutation type is not supported")
        sys.exit(1)


def plot_portrait(sample_dict, percentage=False):
    '''
    Given a dictionary of dataframes, a sample portrait is drawn.
    '''
    with tempfile.TemporaryDirectory() as tempdir:
        os.mkdir(f'{tempdir}/output')
        os.mkdir(f'{tempdir}/output/SBS')
        os.mkdir(f'{tempdir}/output/DBS')
        os.mkdir(f'{tempdir}/output/ID')
        df2csv(df=sample_dict['96'], fname=f'{tempdir}/output/SBS/tmp.SBS96.all')
        df2csv(df=sample_dict['6'], fname=f'{tempdir}/output/SBS/tmp.SBS6.all')
        df2csv(df=sample_dict['24'], fname=f'{tempdir}/output/SBS/tmp.SBS24.all')
        df2csv(df=sample_dict['384'], fname=f'{tempdir}/output/SBS/tmp.SBS384.all')
        df2csv(df=sample_dict['1536'], fname=f'{tempdir}/output/SBS/tmp.SBS1536.all')
        df2csv(df=sample_dict['DBS186'], fname=f'{tempdir}/output/DBS/tmp.DBS186.all')
        df2csv(df=sample_dict['DINUC'], fname=f'{tempdir}/output/DBS/tmp.DBS78.all')
        df2csv(df=sample_dict['ID'], fname=f'{tempdir}/output/ID/tmp.ID83.all')
        df2csv(df=sample_dict['ID415'], fname=f'{tempdir}/output/ID/tmp.ID415.all')
        df2csv(df=sample_dict['ID28'], fname=f'{tempdir}/output/ID/tmp.ID28.all')

        sigPlt_portrait.samplePortrait(f'{tempdir}/', f'{tempdir}/', 'tmp', percentage)
        plot_pdf = convert_from_path(f'{tempdir}/sample_portrait_tmp.pdf', dpi=200)

    return plot_pdf[0]


def _plot_png(sample_name, sig_names, weights, nonzero_exposures, mtype, sample_plot, basis_plot, reconstruction_plot, statistics):
    '''
    result[0].convert('RGB').save("test.pdf")
    result[0].save("test.png")
    '''

    original_png = Image.open(sample_plot[sample_name])
    reconstructed_png = Image.open(reconstruction_plot[sample_name])
    bracket_png = Image.open(paths+"/src/Accolade_fermante.png")
    bases_pngs = [Image.open(v) for k, v in basis_plot.items()]
    bases_order = list(np.argsort(nonzero_exposures))[::-1]
    bases_pngs = [bases_pngs[i] for i in bases_order] 
    bases_combined = Image.fromarray(np.vstack((np.asarray(i) for i in bases_pngs)))
    bases_labels = [weights[i] for i in bases_order]
    summary_string_1 = f"Cosine Similarity: {statistics.iloc[0,0]:.2f}  Correlation: {1 - statistics.iloc[0,2]:.2f}  KL Divergence: {statistics.iloc[0,8]:.2f}"
    summary_string_2 = f"L1 Error %: {statistics.iloc[0,5]:.2f}  L2 Error %: {statistics.iloc[0,7]:.2f}"

    default_width = original_png.size[0]
    default_height = original_png.size[1]

    width = round(default_width * 1.1 + bases_combined.size[0])
    height_left = round(default_height * 2.5)
    height_right = bases_combined.size[1]
    height = max(height_left, height_right)

    if height_right > height_left:
        height_mid = round(height_right/2)
        original_coord = (0, height_mid - default_height)
        reconstructed_coord = (0, height_mid)
        bases_coord = (width - bases_combined.size[0], 0)
        bracket_coord = (default_width, 0)
        bracket_size = (round(default_width * 0.1), height)
    else:
        height_mid = default_height
        original_coord = (0, 0)
        reconstructed_coord = (0, height_mid)
        bases_coord = ((width - bases_combined.size[0], height_mid - round(bases_combined.size[1]/2)))
        bracket_coord = (default_width, 0)
        bracket_size = (round(default_width * 0.1), default_height * 2)

    final_img = Image.new(mode='RGBA', size=(width, height), color=(255,255,255,255))
    draw = ImageDraw.Draw(final_img)
    final_img.paste(original_png, box=original_coord)
    final_img.paste(reconstructed_png, box=reconstructed_coord)
    final_img.paste(bracket_png.resize(bracket_size), box=bracket_coord)
    final_img.paste(bases_combined, box=bases_coord)

    txt_shift_width = round(default_width * 0.80)
    txt_shift_height = round(default_height * 0.2)
    font = ImageFont.truetype("/home/users/kimin/projects/03_Mucosal_Melanoma/scripts/arialbd.ttf", 100)
    draw.text((original_coord[0] + txt_shift_width, original_coord[1] + txt_shift_height), "Original", font=font, fill=(0,0,0,255))
    draw.text((reconstructed_coord[0] + txt_shift_width, reconstructed_coord[1] + txt_shift_height), "Reconstructed", font=font, fill=(0,0,0,255))
    bases_n = len(bases_labels)
    offset_y = round(bases_combined.size[1] / bases_n)
    txt_shift_width = round(bases_combined.size[0] * 0.85)
    for i, text in enumerate(bases_labels):
        draw.text((bases_coord[0] + txt_shift_width, bases_coord[1] + txt_shift_height + offset_y * i), text, font=font, fill=(0,0,0,255))

    font = ImageFont.truetype("/home/users/kimin/projects/03_Mucosal_Melanoma/scripts/arialbd.ttf", 120)
    draw.text((0, round(height*0.85)), summary_string_1, font=font, fill=(0,0,0,255))
    draw.text((0, round(height*0.9)), summary_string_2, font=font, fill=(0,0,0,255))
    
    return final_img


def _plot_decomposition(sample_DF, sample_name, sigDatabases_DF, sig_names, weights, nonzero_exposures, mtype):
    sigAbsoluteCount_DF = sigDatabases_DF[sig_names] * nonzero_exposures
    sigDatabases_DF = pd.concat([sigDatabases_DF.iloc[:,0], sigAbsoluteCount_DF], axis=1)
    sigReconstructed_DF = pd.DataFrame({sample_name:np.sum(sigAbsoluteCount_DF, axis=1)})
    sigReconstructed_DF = pd.concat([sigDatabases_DF.iloc[:,0], sigReconstructed_DF], axis=1)

    with tempfile.TemporaryDirectory() as tmpdirname:
        if mtype == "96":
            sample_plot = sigPlt.plotSBS(sample_DF, tmpdirname, 'test', mtype, False)
            basis_plot = sigPlt.plotSBS(sigDatabases_DF, tmpdirname, 'test', mtype, True)
            reconstruction_plot = sigPlt.plotSBS(sigReconstructed_DF, tmpdirname, "test", mtype, False)
        elif mtype == "78":
            sample_plot = sigPlt.plotDBS(sample_DF, tmpdirname, 'test', mtype, False)
            basis_plot = sigPlt.plotDBS(sigDatabases_DF, tmpdirname, 'test', mtype, True)
            reconstruction_plot = sigPlt.plotDBS(sigReconstructed_DF, tmpdirname, "test", mtype, False)
        elif mtype == "83":
            sample_plot = sigPlt.plotID(sample_DF, tmpdirname, 'test', mtype, False)
            basis_plot = sigPlt.plotID(sigDatabases_DF, tmpdirname, 'test', mtype, True)
            reconstruction_plot = sigPlt.plotID(sigReconstructed_DF, tmpdirname, "test", mtype, False)
        else:
            print("mtype error")
            sys.exit(1)

    statistics = sp.calculate_similarities(sample_DF, sample_name, sigReconstructed_DF[sample_name]/np.sum(nonzero_exposures))

    img = _plot_png(sample_name, sig_names, weights, nonzero_exposures, mtype, sample_plot, basis_plot, reconstruction_plot, statistics)

    return img


def sample_decomposition(sample, signatures=None, mtype="", genome_build="GRCh37", add_penalty=0.05, remove_penalty=0.01, mutation_context=None, connected_sigs=True, make_decomposition_plots=True, originalProcessAvg=None):
    '''
    Wrapper function edited from SigProfilerExtractor.subroutines.signature_decomposition
    Takes one dataframe of original mutation count, and optionally takes a dataframe of custom signatures other than cosmic (signatures, default None).
    Returns a tuple containing PIL img, decomposed signature df, exposure df
    '''

    if sample.shape[0]==96:
        mtype = '96'
        if genome_build=="GRCh37":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="GRCh37")
            sigDatabase,_,_,_ = sub.read_csv(sigDatabase)
        elif genome_build=="GRCh38":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="GRCh38")
            sigDatabase,_,_,_ = sub.read_csv(sigDatabase) 
        elif genome_build=="mm9":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="mm9")
            sigDatabase,_,_,_ = sub.read_csv(sigDatabase)  
        elif genome_build=="mm10":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="mm10")
            sigDatabase,_,_,_ = sub.read_csv(sigDatabase) 
            
        signames = sigDatabase.columns 


    elif sample.shape[0]==78:
        mtype = '78'
        if genome_build=="GRCh37":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="GRCh37", index_col=0)
        elif genome_build=="GRCh38":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="GRCh38", index_col=0)
        elif genome_build=="mm9":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="mm9",index_col=0)
        elif genome_build=="mm10":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="mm10",index_col=0)
            
        signames = sigDatabase.columns
        connected_sigs=False


    elif sample.shape[0]==83:
        mtype = '83'
        sigDatabase = pd.read_csv(paths+"/data/sigProfiler_ID_signatures.csv", sep=",",index_col=0)
        signames = sigDatabase.columns
        connected_sigs=False
        

    elif sample.shape[0]==48:
        mtype = '48'
        sigDatabase = pd.read_csv(paths+"/data/CNV_signatures.txt", sep="\t",index_col=0)
        signames = sigDatabase.columns
        connected_sigs=False
        print(list(signames))
        print(list(sigDatabase.index))


    else:
        mtype = str(sample.shape[0])
        print(f"Customized signatures' mtype is {mtype}")
        sigDatabase = pd.DataFrame(signatures)
        sigDatabase.columns=sigDatabase.columns.astype(str)
        sigDatabase.index=sigDatabase.index.astype(str)
        signames=sigDatabase.columns
        connected_sigs=False


    if not np.all(np.equal(sample.index, sigDatabase.index)):
        print("sample index and signature index are different")
        sys.exit(1)
    

    sigDatabase_array = np.array(sigDatabase)
    sample_array = np.array(sample).ravel()


    ## If mutation context is 96, the fit process is little different 
    if mtype=='96':
        if genome_build=="mm9" or genome_build=="mm10":
            check_rule_negatives = [1,16]
            check_rule_penalty=1.50
        else:
            check_rule_negatives = []
            check_rule_penalty=1.0

        with tempfile.TemporaryDirectory() as tmpdirname:
            _, exposures,L2dist,similarity, kldiv, correlation, cosine_similarity_with_four_signatures = ss.add_remove_signatures(sigDatabase_array, 
                                                                                             sample_array, 
                                                                                             metric="l2", 
                                                                                             solver="nnls", 
                                                                                             background_sigs = [0,4], 
                                                                                             permanent_sigs = [0,4], 
                                                                                             candidate_sigs="all", 
                                                                                             allsigids = signames, 
                                                                                             add_penalty = add_penalty, 
                                                                                             remove_penalty = remove_penalty,
                                                                                             check_rule_negatives = check_rule_negatives, 
                                                                                             checkrule_penalty = check_rule_penalty, 
                                                                                             directory = tmpdirname+"/log.txt", 
                                                                                             connected_sigs=connected_sigs,
                                                                                             verbose=False)
            with open(tmpdirname+"/log.txt", 'r') as log:
                print(log.read())


    ## If mutation context is not 96
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            _, exposures,L2dist,similarity, kldiv, correlation, cosine_similarity_with_four_signatures = ss.add_remove_signatures(sigDatabase_array, 
                                                                                             sample_array, 
                                                                                             metric="l2", 
                                                                                             solver="nnls", 
                                                                                             background_sigs = [], 
                                                                                             permanent_sigs = [], 
                                                                                             candidate_sigs="all",
                                                                                             allsigids = signames,  
                                                                                             add_penalty = add_penalty, 
                                                                                             remove_penalty = remove_penalty,
                                                                                             check_rule_negatives = [], 
                                                                                             checkrule_penalty = [], 
                                                                                             directory = tmpdirname+"/log.txt", 
                                                                                             connected_sigs=connected_sigs,
                                                                                             verbose=False)
            with open(tmpdirname+"/log.txt", 'r') as log:
                print(log.read())


    # calculate the L1 Error %
    L1dist = np.linalg.norm(sample_array-np.dot(sigDatabase,exposures) , ord=1)/np.linalg.norm(sample_array, ord=1)
    exposure_percentages = exposures[np.nonzero(exposures)]/np.sum(exposures[np.nonzero(exposures)])*100

    count =0
    listofinformation = list("0"*len(np.nonzero(exposures)[0])*3)
    decomposed_signatures = []
    contribution_percentages = []
    
    for j in np.nonzero(exposures)[0]:
        listofinformation[count*3] = signames[j]
        decomposed_signatures.append(signames[j])
        listofinformation[count*3+1] = round(exposure_percentages[count],2)
        contribution_percentages.append(round(exposure_percentages[count],2))
        listofinformation[count*3+2]="%"
        count+=1

    weights=[]
    basis_names=[]
    nonzero_exposures=exposures[np.nonzero(exposures)]
    for info in range(0, len(listofinformation), 3):
        sigName=listofinformation[info]
        basis_names.append(sigName)
        sigWeigt=str(listofinformation[info+1])+"%"
        weights.append(sigWeigt)



    ## Plotting preparation
    sample_DF = sample.reset_index()
    sample_name = sample.columns[0]

    sigDatabases_DF = sigDatabase.reset_index()
    sigDatabases_name = sigDatabases_DF.columns[0]

    basis_cols = basis_names.copy()
    basis_cols.insert(0,sigDatabases_name)

    img = _plot_decomposition(sample_DF, sample_name, sigDatabases_DF[basis_cols], basis_names, weights, nonzero_exposures, mtype)
    # resize image
    w,h = img.size
    img = img.resize((round(w/6), round(h/6)))

    sig_summary = pd.DataFrame({'Signature': basis_names, 'exposure': nonzero_exposures, 'prop': weights}).sort_values(by='exposure', ascending=False)

    return img, sigDatabase[basis_names], sig_summary
    



def sample_decomposition_old(sample, signatures=None, mtype="", genome_build="GRCh37", add_penalty=0.05, remove_penalty=0.01, mutation_context=None, connected_sigs=True, make_decomposition_plots=True, originalProcessAvg=None):
    '''
    Wrapper function edited from SigProfilerExtractor.subroutines.signature_decomposition
    Takes one dataframe of original mutation count, and optionally takes a dataframe of custom signatures other than cosmic (signatures, default None).
    Returns a dictionary containing PIL img, decomposed signature ids, decomposed signature matrix, exposure array etc.
    The process is inefficient as the img is created by pdf, saved in a temp folder, and loaded back again by PIL module.
    '''
    paths = cosmic.__path__[0]
    merger = PdfFileMerger()

    if sample.shape[0]==96:
        mtype = '96'
        if genome_build=="GRCh37":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="GRCh37")
            sigDatabase,_,_,_ = sub.read_csv(sigDatabase)
        elif genome_build=="GRCh38":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="GRCh38")
            sigDatabase,_,_,_ = sub.read_csv(sigDatabase) 
        elif genome_build=="mm9":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="mm9")
            sigDatabase,_,_,_ = sub.read_csv(sigDatabase)  
        elif genome_build=="mm10":
            sigDatabase = pd.read_excel(paths+"/data/SBS_signatures_genome_builds.xlsx", sheet_name="mm10")
            sigDatabase,_,_,_ = sub.read_csv(sigDatabase) 
            
        signames = sigDatabase.columns 


    elif sample.shape[0]==78:
        mtype = '78'
        if genome_build=="GRCh37":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="GRCh37", index_col=0)
        elif genome_build=="GRCh38":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="GRCh38", index_col=0)
        elif genome_build=="mm9":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="mm9",index_col=0)
        elif genome_build=="mm10":
            sigDatabase = pd.read_excel(paths+"/data/DBS_signatures_genome_builds.xlsx", sheet_name="mm10",index_col=0)
            
        signames = sigDatabase.columns
        connected_sigs=False


    elif sample.shape[0]==83:
        mtype = '83'
        sigDatabase = pd.read_csv(paths+"/data/sigProfiler_ID_signatures.csv", sep=",",index_col=0)
        signames = sigDatabase.columns
        connected_sigs=False
        

    elif sample.shape[0]==48:
        mtype = '48'
        sigDatabase = pd.read_csv(paths+"/data/CNV_signatures.txt", sep="\t",index_col=0)
        signames = sigDatabase.columns
        connected_sigs=False
        print(list(signames))
        print(list(sigDatabase.index))


    else:
        mtype = str(sample.shape[0])
        print(f"Customized signatures' mtype is {mtype}")
        sigDatabase = pd.DataFrame(signatures)
        sigDatabase.columns=sigDatabase.columns.astype(str)
        sigDatabase.index=sigDatabase.index.astype(str)
        signames=sigDatabase.columns
        connected_sigs=False


    if not np.all(np.equal(sample.index, sigDatabase.index)):
        print("sample index and signature index are different")
        sys.exit(1)
    

    sigDatabase_array = np.array(sigDatabase)
    sample_array = np.array(sample).ravel()


    ## If mutation context is 96, the fit process is little different 
    if mtype=='96':
        if genome_build=="mm9" or genome_build=="mm10":
            check_rule_negatives = [1,16]
            check_rule_penalty=1.50
        else:
            check_rule_negatives = []
            check_rule_penalty=1.0

        with tempfile.TemporaryDirectory() as tmpdirname:
            _, exposures,L2dist,similarity, kldiv, correlation, cosine_similarity_with_four_signatures = ss.add_remove_signatures(sigDatabase_array, 
                                                                                             sample_array, 
                                                                                             metric="l2", 
                                                                                             solver="nnls", 
                                                                                             background_sigs = [0,4], 
                                                                                             permanent_sigs = [0,4], 
                                                                                             candidate_sigs="all", 
                                                                                             allsigids = signames, 
                                                                                             add_penalty = add_penalty, 
                                                                                             remove_penalty = remove_penalty,
                                                                                             check_rule_negatives = check_rule_negatives, 
                                                                                             checkrule_penalty = check_rule_penalty, 
                                                                                             directory = tmpdirname+"/log.txt", 
                                                                                             connected_sigs=connected_sigs,
                                                                                             verbose=False)
            with open(tmpdirname+"/log.txt", 'r') as log:
                print(log.read())


    ## If mutation context is not 96
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            _, exposures,L2dist,similarity, kldiv, correlation, cosine_similarity_with_four_signatures = ss.add_remove_signatures(sigDatabase_array, 
                                                                                             sample_array, 
                                                                                             metric="l2", 
                                                                                             solver="nnls", 
                                                                                             background_sigs = [], 
                                                                                             permanent_sigs = [], 
                                                                                             candidate_sigs="all",
                                                                                             allsigids = signames,  
                                                                                             add_penalty = add_penalty, 
                                                                                             remove_penalty = remove_penalty,
                                                                                             check_rule_negatives = [], 
                                                                                             checkrule_penalty = [], 
                                                                                             directory = tmpdirname+"/log.txt", 
                                                                                             connected_sigs=connected_sigs,
                                                                                             verbose=False)
            with open(tmpdirname+"/log.txt", 'r') as log:
                print(log.read())


    # calculate the L1 Error %
    L1dist = np.linalg.norm(sample_array-np.dot(sigDatabase,exposures) , ord=1)/np.linalg.norm(sample_array, ord=1)
    exposure_percentages = exposures[np.nonzero(exposures)]/np.sum(exposures[np.nonzero(exposures)])*100

    count =0
    listofinformation = list("0"*len(np.nonzero(exposures)[0])*3)
    decomposed_signatures = []
    contribution_percentages = []
    
    for j in np.nonzero(exposures)[0]:
        listofinformation[count*3] = signames[j]
        decomposed_signatures.append(signames[j])
        listofinformation[count*3+1] = round(exposure_percentages[count],2)
        contribution_percentages.append(round(exposure_percentages[count],2))
        listofinformation[count*3+2]="%"
        count+=1

    weights=[]
    basis_names=[]
    nonzero_exposures=exposures[np.nonzero(exposures)]
    for info in range(0, len(listofinformation), 3):
        sigName=listofinformation[info]
        basis_names.append(sigName)
        sigWeigt=str(listofinformation[info+1])+"%"
        weights.append(sigWeigt)



    ## Plotting preparation
    sample_DF = sample.reset_index()
    sample_name = sample.columns[0]

    sigDatabases_DF = sigDatabase.reset_index()
    sigDatabases_name = sigDatabases_DF.columns[0]

    basis_cols = basis_names.copy()
    basis_cols.insert(0,sigDatabases_name)

     
    with tempfile.TemporaryDirectory() as tmpdirname:
        byte_plot = sp.run_PlotDecomposition(sample_DF, sample_name, sigDatabases_DF[basis_cols], basis_names, weights, exposure_percentages/100, tmpdirname, "tmp", mtype)
        merger.append(byte_plot)
        merger.write(tmpdirname+"/test_Decomposition_Plots.pdf")
        img = convert_from_path(tmpdirname+"/test_Decomposition_Plots.pdf", dpi=200)[0] ## Increase dpi for better image resolution. 300 is okay for presentation I think. 

    
    print("Decompositon Plot made for {}\n".format(sample_name))



    ## Process data for return
    different_signatures = np.unique(np.nonzero(exposures))
    different_signatures=different_signatures.astype(int)
    if mtype == "96" or mtype=="288" or mtype=="1536":
        different_signatures = list(set().union(different_signatures, [0,4]))
        different_signatures.sort()    
      
    
    #get the name of the signatures
    try:
        detected_signatures = signames[different_signatures]
        globalsigmats= sigDatabase.loc[:,list(detected_signatures)]
    except:
        detected_signatures=[None]
        globalsigmats=None
           
    #only for SBS96
    if mtype == "96" or mtype=="288" or mtype=="1536":        
        background_sigs = sub.get_indeces(list(detected_signatures), ['SBS1', 'SBS5'])
        # add connected signatures   
        different_signatures = ss.add_connected_sigs(different_signatures, list(signames))
    #for other contexts
    else:
        background_sigs = []


    return {"img": img, "globalsigids": list(detected_signatures), "globalsigs":globalsigmats,  
            "background_sigs": background_sigs, "contribution_percentages": contribution_percentages} 



def single_decomposition(vcf_path, output=None, ref='GRCh37', sig_database="default", exome=False):
    '''
    Wrapper function edited from SigProfilerSingleSample
    The process is little bit different from SigProfilerExtract.decompose. The difference is from the argument check_rules. However this check_rules does not apply to ID or DBS decomposition.
    For 96, it is okay to use this wrapper function, but for other mutation context, use sample_decomposition. 
    
    vcf_path: The original function can take df and vcf_path as as input, but using df here set check_rules to false. If check_rules is false, it's better to use sample_decomposition.
    output: Output folder path just in case you want to save the outcome. If none, temp folder will be created.
    sig_database: Custom df, default "default"
    exome: Default False
    '''

    if not output:
        tmpdir = tempfile.TemporaryDirectory()
        output = tmpdir.name
    else:
        output = output.strip("/")
        subprocess.run(f"mkdir -p {output}")
    print("output folder is ", output)


    run = subprocess.run(f'cp {vcf_path} {output} ', shell=True, capture_output=True)
    # print(run.stdout)
    # print(run.stderr)
    if vcf_path.endswith('gz'):
        run = subprocess.run(f'gzip -d {output}/{os.path.basename(vcf_path)}', shell=True, capture_output=True)
        # print(run.stdout)
        # print(run.stderr)
        input = f'{output}/{os.path.basename(vcf_path).strip(".gz")}'
    else:
        input = f'{output}/{os.path.basename(vcf_path)}'
    print("input file is ", input)

    spss.single_sample(output, output, ref=ref, sig_database=sig_database, check_rules=True, exome=exome)

    sample = pd.read_csv(f'{output}/output/SBS/{os.path.basename(output)}.SBS96.all', index_col=0, sep="\t")
    sigDatabase = pd.read_csv(f'{output}/Signatures.txt', index_col=0, sep="\t")
    nonzero_exposures = np.genfromtxt(f'{output}/Sig_activities.txt', delimiter='\t', skip_header=1, usecols=1, dtype=int)
    exposure_percentages = nonzero_exposures/np.sum(nonzero_exposures)*100
    weights = [f"{x:.2f}%" for x in exposure_percentages]
    probabilities = pd.read_csv(f'{output}/Mutation_Probabilities.txt', index_col=1, sep="\t")
    probabilities = probabilities.drop(probabilities.columns[0], axis=1)

    sample_DF = sample.reset_index()
    sample_name = sample.columns[0]
    mtype = str(sample.shape[0])

    sigDatabases_DF = sigDatabase.reset_index()
    basis_names = sigDatabase.columns

    img = _plot_decomposition(sample_DF, sample_name, sigDatabases_DF, basis_names, weights, nonzero_exposures, mtype)

    sig_summary = pd.DataFrame({'Signature': basis_names, 'exposure': nonzero_exposures, 'prop': weights}).sort_values(by='exposure', ascending=False)

    return img, sigDatabase, sig_summary, probabilities




## Sig profiler 에서 가져온 것들
## 안쓰더라도 나중에 유용할까봐 가져옴
# Provides a chromosome conversion from NCBI notation


    # letters = list(string.ascii_uppercase)
    # letters.extend([i+b for i in letters for b in letters])
    # letters = letters[0:signatures.shape[1]]



def perm(n, seq):
    '''
    From Sig profiler

    Generates a list of all available permutations of n-mers.
    Parameters:
               n  -> length of the desired permutation string
             seq  -> list of all possible string valuesvcf_pathvcf_path
    Returns:
          permus  -> list of all available permutations
    '''
    permus = []
    for p in itertools.product(seq, repeat=n):
        permus.append("".join(p))
    return(permus)




bases = ['A','C','G','T']
tsb = ['T','U','N','B']

ncbi_chrom = {'NC_000067.6':'1', 'NC_000068.7':'2', 'NC_000069.6':'3', 'NC_000070.6':'4', 
              'NC_000071.6':'5', 'NC_000072.6':'6', 'NC_000073.6':'7', 'NC_000074.6':'8',
              'NC_000075.6':'9', 'NC_000076.6':'10', 'NC_000077.6':'11', 'NC_000078.6':'12',
              'NC_000079.6':'13', 'NC_000080.6':'14', 'NC_000081.6':'15', 'NC_000082.6':'16', 
              'NC_000083.6':'17', 'NC_000084.6':'18', 'NC_000085.6':'19', 'NC_000086.7':'X', 
              'NC_000087.7':'Y', '82503188|ref|NC_007605.1|':'gi_82503188_ref_NC_007605'}

mutation_types = ['CC>AA','CC>AG','CC>AT','CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT',
              'CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG','CT>TA','CT>TC','CT>TG',
              'TC>AA','TC>AG','TC>AT','TC>CA','TC>CG','TC>CT','TC>GA','TC>GG','TC>GT',
              'TT>AA','TT>AC','TT>AG','TT>CA','TT>CC','TT>CG','TT>GA','TT>GC','TT>GG']

mutation_types_non_tsb = ['AC>CA','AC>CG','AC>CT','AC>GA','AC>GG','AC>GT','AC>TA','AC>TG','AC>TT',
              'AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA',
              'CG>AT','CG>GC','CG>GT','CG>TA','CG>TC','CG>TT',
              'GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA',
              'TA>AT','TA>CG','TA>CT','TA>GC','TA>GG','TA>GT',
              'TG>AA','TG>AC','TG>AT','TG>CA','TG>CC','TG>CT','TG>GA','TG>GC','TG>GT']

indels_seq_types = [ # Single-sequences
                'C', 'T',

                # Di-sequences
                'AC','AT','CA','CC','CG','CT','GC','TA','TC','TT',

                # Tri-sequences
                'ACC', 'ACT', 'ATC', 'ATT', 'CAC', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGC', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 
                'GCC', 'GCT', 'GTC', 'GTT', 'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT',

                # Tetra-sequences
                'AACC', 'AACT', 'AATC', 'AATT', 'ACAC', 'ACAT', 'ACCA', 'ACCC', 'ACCG', 'ACCT', 'ACGC', 'ACGT', 'ACTA', 'ACTC', 'ACTG', 'ACTT', 'AGCC', 'AGCT', 'AGTC', 
                 'AGTT', 'ATAC', 'ATAT', 'ATCA', 'ATCC', 'ATCG', 'ATCT', 'ATGC', 'ATGT', 'ATTA', 'ATTC', 'ATTG', 'ATTT', 'CAAC', 'CAAT', 'CACA', 'CACC', 'CACG', 'CACT', 
                 'CAGC', 'CAGT', 'CATA', 'CATC', 'CATG', 'CATT', 'CCAA', 'CCAC', 'CCAG', 'CCAT', 'CCCA', 'CCCC', 'CCCG', 'CCCT', 'CCGA', 'CCGC', 'CCGG', 'CCGT', 'CCTA', 
                 'CCTC', 'CCTG', 'CCTT', 'CGAC', 'CGAT', 'CGCA', 'CGCC', 'CGCG', 'CGCT', 'CGGC', 'CGTA', 'CGTC', 'CGTG', 'CGTT', 'CTAA', 'CTAC', 'CTAG', 'CTAT', 'CTCA', 
                 'CTCC', 'CTCG', 'CTCT', 'CTGA', 'CTGC', 'CTGG', 'CTGT', 'CTTA', 'CTTC', 'CTTG', 'CTTT', 'GACC', 'GATC', 'GCAC', 'GCCA', 'GCCC', 'GCCG', 'GCCT', 'GCGC', 
                 'GCTA', 'GCTC', 'GCTG', 'GCTT', 'GGCC', 'GGTC', 'GTAC', 'GTCA', 'GTCC', 'GTCG', 'GTCT', 'GTGC', 'GTTA', 'GTTC', 'GTTG', 'GTTT', 'TAAC', 'TACA', 'TACC', 
                 'TACG', 'TACT', 'TAGC', 'TATA', 'TATC', 'TATG', 'TATT', 'TCAA', 'TCAC', 'TCAG', 'TCAT', 'TCCA', 'TCCC', 'TCCG', 'TCCT', 'TCGA', 'TCGC', 'TCGG', 'TCGT', 
                 'TCTA', 'TCTC', 'TCTG', 'TCTT', 'TGAC', 'TGCA', 'TGCC', 'TGCG', 'TGCT', 'TGTA', 'TGTC', 'TGTG', 'TGTT', 'TTAA', 'TTAC', 'TTAG', 'TTAT', 'TTCA', 'TTCC', 
                 'TTCG', 'TTCT', 'TTGA', 'TTGC', 'TTGG', 'TTGT', 'TTTA', 'TTTC', 'TTTG', 'TTTT',

                # Penta-sequences
                'AACCC', 'AACCT', 'AACTC', 'AACTT', 'AATCC', 'AATCT', 'AATTC', 'AATTT', 'ACACC', 'ACACT', 'ACATC', 'ACATT', 'ACCAC', 'ACCAT', 'ACCCA', 'ACCCC', 'ACCCG', 
                'ACCCT', 'ACCGC', 'ACCGT', 'ACCTA', 'ACCTC', 'ACCTG', 'ACCTT', 'ACGCC', 'ACGCT', 'ACGTC', 'ACGTT', 'ACTAC', 'ACTAT', 'ACTCA', 'ACTCC', 'ACTCG', 'ACTCT', 
                'ACTGC', 'ACTGT', 'ACTTA', 'ACTTC', 'ACTTG', 'ACTTT', 'AGCCC', 'AGCCT', 'AGCTC', 'AGCTT', 'AGTCC', 'AGTCT', 'AGTTC', 'AGTTT', 'ATACC', 'ATACT', 'ATATC', 
                'ATATT', 'ATCAC', 'ATCAT', 'ATCCA', 'ATCCC', 'ATCCG', 'ATCCT', 'ATCGC', 'ATCGT', 'ATCTA', 'ATCTC', 'ATCTG', 'ATCTT', 'ATGCC', 'ATGCT', 'ATGTC', 'ATGTT', 
                'ATTAC', 'ATTAT', 'ATTCA', 'ATTCC', 'ATTCG', 'ATTCT', 'ATTGC', 'ATTGT', 'ATTTA', 'ATTTC', 'ATTTG', 'ATTTT', 'CAACC', 'CAACT', 'CAATC', 'CAATT', 'CACAC', 
                'CACAT', 'CACCA', 'CACCC', 'CACCG', 'CACCT', 'CACGC', 'CACGT', 'CACTA', 'CACTC', 'CACTG', 'CACTT', 'CAGCC', 'CAGCT', 'CAGTC', 'CAGTT', 'CATAC', 'CATAT', 
                'CATCA', 'CATCC', 'CATCG', 'CATCT', 'CATGC', 'CATGT', 'CATTA', 'CATTC', 'CATTG', 'CATTT', 'CCAAC', 'CCAAT', 'CCACA', 'CCACC', 'CCACG', 'CCACT', 'CCAGC', 
                'CCAGT', 'CCATA', 'CCATC', 'CCATG', 'CCATT', 'CCCAA', 'CCCAC', 'CCCAG', 'CCCAT', 'CCCCA', 'CCCCC', 'CCCCG', 'CCCCT', 'CCCGA', 'CCCGC', 'CCCGG', 'CCCGT', 
                'CCCTA', 'CCCTC', 'CCCTG', 'CCCTT', 'CCGAC', 'CCGAT', 'CCGCA', 'CCGCC', 'CCGCG', 'CCGCT', 'CCGGC', 'CCGGT', 'CCGTA', 'CCGTC', 'CCGTG', 'CCGTT', 'CCTAA', 
                'CCTAC', 'CCTAG', 'CCTAT', 'CCTCA', 'CCTCC', 'CCTCG', 'CCTCT', 'CCTGA', 'CCTGC', 'CCTGG', 'CCTGT', 'CCTTA', 'CCTTC', 'CCTTG', 'CCTTT', 'CGACC', 'CGACT', 
                'CGATC', 'CGATT', 'CGCAC', 'CGCAT', 'CGCCA', 'CGCCC', 'CGCCG', 'CGCCT', 'CGCGC', 'CGCGT', 'CGCTA', 'CGCTC', 'CGCTG', 'CGCTT', 'CGGCC', 'CGGCT', 'CGGTC', 
                'CGGTT', 'CGTAC', 'CGTAT', 'CGTCA', 'CGTCC', 'CGTCG', 'CGTCT', 'CGTGC', 'CGTGT', 'CGTTA', 'CGTTC', 'CGTTG', 'CGTTT', 'CTAAC', 'CTAAT', 'CTACA', 'CTACC', 
                'CTACG', 'CTACT', 'CTAGC', 'CTAGT', 'CTATA', 'CTATC', 'CTATG', 'CTATT', 'CTCAA', 'CTCAC', 'CTCAG', 'CTCAT', 'CTCCA', 'CTCCC', 'CTCCG', 'CTCCT', 'CTCGA', 
                'CTCGC', 'CTCGG', 'CTCGT', 'CTCTA', 'CTCTC', 'CTCTG', 'CTCTT', 'CTGAC', 'CTGAT', 'CTGCA', 'CTGCC', 'CTGCG', 'CTGCT', 'CTGGC', 'CTGGT', 'CTGTA', 'CTGTC', 
                'CTGTG', 'CTGTT', 'CTTAA', 'CTTAC', 'CTTAG', 'CTTAT', 'CTTCA', 'CTTCC', 'CTTCG', 'CTTCT', 'CTTGA', 'CTTGC', 'CTTGG', 'CTTGT', 'CTTTA', 'CTTTC', 'CTTTG', 
                'CTTTT', 'GACCC', 'GACCT', 'GACTC', 'GACTT', 'GATCC', 'GATCT', 'GATTC', 'GATTT', 'GCACC', 'GCACT', 'GCATC', 'GCATT', 'GCCAC', 'GCCAT', 'GCCCA', 'GCCCC', 
                'GCCCG', 'GCCCT', 'GCCGC', 'GCCGT', 'GCCTA', 'GCCTC', 'GCCTG', 'GCCTT', 'GCGCC', 'GCGCT', 'GCGTC', 'GCGTT', 'GCTAC', 'GCTAT', 'GCTCA', 'GCTCC', 'GCTCG', 
                'GCTCT', 'GCTGC', 'GCTGT', 'GCTTA', 'GCTTC', 'GCTTG', 'GCTTT', 'GGCCC', 'GGCCT', 'GGCTC', 'GGCTT', 'GGTCC', 'GGTCT', 'GGTTC', 'GGTTT', 'GTACC', 'GTACT', 
                'GTATC', 'GTATT', 'GTCAC', 'GTCAT', 'GTCCA', 'GTCCC', 'GTCCG', 'GTCCT', 'GTCGC', 'GTCGT', 'GTCTA', 'GTCTC', 'GTCTG', 'GTCTT', 'GTGCC', 'GTGCT', 'GTGTC', 
                'GTGTT', 'GTTAC', 'GTTAT', 'GTTCA', 'GTTCC', 'GTTCG', 'GTTCT', 'GTTGC', 'GTTGT', 'GTTTA', 'GTTTC', 'GTTTG', 'GTTTT', 'TAACC', 'TAACT', 'TAATC', 'TAATT', 
                'TACAC', 'TACAT', 'TACCA', 'TACCC', 'TACCG', 'TACCT', 'TACGC', 'TACGT', 'TACTA', 'TACTC', 'TACTG', 'TACTT', 'TAGCC', 'TAGCT', 'TAGTC', 'TAGTT', 'TATAC', 
                'TATAT', 'TATCA', 'TATCC', 'TATCG', 'TATCT', 'TATGC', 'TATGT', 'TATTA', 'TATTC', 'TATTG', 'TATTT', 'TCAAC', 'TCAAT', 'TCACA', 'TCACC', 'TCACG', 'TCACT', 
                'TCAGC', 'TCAGT', 'TCATA', 'TCATC', 'TCATG', 'TCATT', 'TCCAA', 'TCCAC', 'TCCAG', 'TCCAT', 'TCCCA', 'TCCCC', 'TCCCG', 'TCCCT', 'TCCGA', 'TCCGC', 'TCCGG', 
                'TCCGT', 'TCCTA', 'TCCTC', 'TCCTG', 'TCCTT', 'TCGAC', 'TCGAT', 'TCGCA', 'TCGCC', 'TCGCG', 'TCGCT', 'TCGGC', 'TCGGT', 'TCGTA', 'TCGTC', 'TCGTG', 'TCGTT', 
                'TCTAA', 'TCTAC', 'TCTAG', 'TCTAT', 'TCTCA', 'TCTCC', 'TCTCG', 'TCTCT', 'TCTGA', 'TCTGC', 'TCTGG', 'TCTGT', 'TCTTA', 'TCTTC', 'TCTTG', 'TCTTT', 'TGACC', 
                'TGACT', 'TGATC', 'TGATT', 'TGCAC', 'TGCAT', 'TGCCA', 'TGCCC', 'TGCCG', 'TGCCT', 'TGCGC', 'TGCGT', 'TGCTA', 'TGCTC', 'TGCTG', 'TGCTT', 'TGGCC', 'TGGCT', 
                'TGGTC', 'TGGTT', 'TGTAC', 'TGTAT', 'TGTCA', 'TGTCC', 'TGTCG', 'TGTCT', 'TGTGC', 'TGTGT', 'TGTTA', 'TGTTC', 'TGTTG', 'TGTTT', 'TTAAC', 'TTAAT', 'TTACA', 
                'TTACC', 'TTACG', 'TTACT', 'TTAGC', 'TTAGT', 'TTATA', 'TTATC', 'TTATG', 'TTATT', 'TTCAA', 'TTCAC', 'TTCAG', 'TTCAT', 'TTCCA', 'TTCCC', 'TTCCG', 'TTCCT', 
                'TTCGA', 'TTCGC', 'TTCGG', 'TTCGT', 'TTCTA', 'TTCTC', 'TTCTG', 'TTCTT', 'TTGAC', 'TTGAT', 'TTGCA', 'TTGCC', 'TTGCG', 'TTGCT', 'TTGGC', 'TTGGT', 'TTGTA', 
                'TTGTC', 'TTGTG', 'TTGTT', 'TTTAA', 'TTTAC', 'TTTAG', 'TTTAT', 'TTTCA', 'TTTCC', 'TTTCG', 'TTTCT', 'TTTGA', 'TTTGC', 'TTTGG', 'TTTGT', 'TTTTA', 'TTTTC', 
                'TTTTG', 'TTTTT']

size = 5
mut_types_initial = perm(size, "ACGT")
mut_types = []
for tsbs in tsb:
    for mut in mut_types_initial:
        current_base = mut[int(size/2)]
        if current_base == 'C' or current_base == 'T':
            for base in bases:
                if base != current_base:
                    mut_types.append(tsbs+":"+mut[0:int(size/2)] + "[" + current_base+">"+ base+"]"+mut[int(size/2)+1:])

mutation_types_tsb_context = []
for base in bases:
    for mut in mutation_types:
        for base2 in bases:
            for base3 in tsb:
                mutation_types_tsb_context.append(''.join([base3,":",base,"[",mut,"]",base2]))

for base in bases:
    for mut in mutation_types_non_tsb:
        for base2 in bases:
            mutation_types_tsb_context.append(''.join(['Q:', base, "[", mut, "]", base2]))

