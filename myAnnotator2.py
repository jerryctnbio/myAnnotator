"""This script annotates variants in a VCF file. Besides the depth, read count
directly from the VCF file, it also queries the ExAC database
(http://exac.hms.harvard.edu) for minor allele frequency. The variant effects
are by vep running off a local cache offline. Only the most deleterious effect
is kept in the output for each variant in the VCF file, along with the allele
frequency from the ExAC database. This applies whether the variant matches to
multiple transcripts or the variant has multiple minor alleles, or both. 
See the accompanying README for details. 
"""

__version='R1.0.0.0'


import os, sys
from pathlib import Path

import requests
import time

import argparse
import subprocess


def main():

    """define the command line arguments.

    1. input VCF filename, required.
    2. output tsv file name, optional. Defaults to match input file.
       For example, if input='x.vcf', output='x.annotations.tsv'.
    3. filter flag, defaults to **NOT** filter out variants failing FILTER.
       If a variant failed the FILTER the field would be marked with specific
       filter names. Otherwise, it is marked as 'PASS' or '.'.
    4. flag for variant severity file, defaults to use a hard code file name.
       Users can specify their own one here.
       The default file is "myAnnotator.severity.txt", in the 
                 ***current working directory***.
    5. location of local cache for vep.
    6. genome version, whether it is GRCh37 or GRCh38.
    """  

    title="This script annotates a VCF file with additional information from\n"
    title+="ExAC database and vep.\n"
    title+="The information from the VCF file are read depth, read count\n"
    title+="and frequency for all the alleles.\n"
    title+="Each variant is annotated with the most deleterious effect,\n"
    title+="allele frequency from ExAC, and variant class (e.g. 'snv')\n"
    title+="from vep.\n"
    title+=" If applicable, the gene symbols are also returned.\n"

    argParser=argparse.ArgumentParser(description=title)

    argParser.add_argument("-i", required=True 
                               , help="input VCF file name")
    argParser.add_argument("-o"
                               , help="optional, output tsv file name")

    msg="optional, variant severity level file name"
    argParser.add_argument("-s"
                               , help=msg)

    msg="optional, flag to whether to filter out variants failing FILTER"
    argParser.add_argument("-f", default=False 
                               , action="store_true", help=msg)

    msg="optional, directory of local cache, used for vep'"
    argParser.add_argument("-c" 
                               , help=msg)

    msg='optional, select a human genome builds: GRCh37 or GRCh38'
    argParser.add_argument("-g", choices=['GRCh37', 'GRCh38']
                               , help=msg, default='GRCh37')

    args=argParser.parse_args()

    inFile=args.i
    filter=args.f
    severityFile=args.s
    genome=args.g
    cache=args.c

    severity_dict=setupSeverityTable(severityFile)

    vep_vcf_file=runVep(inFile, genome, cache)

    anno=parseVepVcf(severity_dict, vep_vcf_file, filter)

    # save annotations to a file
    if args.o==None:

        parent=Path(inFile).parent
        stem=Path(inFile).stem
        name=stem+'.annotations.tsv'
        outFile=Path(parent/name)

    else:
        outFile=args.o

    saveListToTabFile(anno, outFile)


def runVep(inFile:str, genome:str, cache:str)->str:

    """Run vep to annotate all the variants. The consequences from vep is
    stored in the output vcf file as the item 'CSQ' in the INFO field.

    Parameters:
    inFile      : str  --- inut file name
    genome      : str  --- whether 'GRCh37' or 'GRCh38'
                           Defaults to 'GRCh37' to match the input vcf file. 
    cache       : str  --- directory for a local cache used by vep.
                           Useful if the cache is not in the default place.

    Returns:
    A string    --- the output vcf file name from vep.

    *** vep has to be in the $Path environment for it to work ***
    """

    field_list=[]
    field_list+=['SYMBOL', 'VARIANT_CLASS', 'Consequence']
    field_list+=['Allele']

    field=",".join(field_list)

    outFile="__vep__temp__.out.txt"

    # completely offline
    cmd=['vep', '--force_overwrite', '--cache', '--symbol', '--offline']

    # specify a genome version
    cmd+=['-assembly', genome]

    if cache !=None:
        cmd+=['--dir_cache', cache]
   
    # output a vcf file 
    cmd+=['--vcf']

    cmd+=['--fields', field, '--variant_class']
    cmd+=['-i', inFile, '-o', outFile]

    print(' '.join(cmd))

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
       
        print(f"\nvep failed. The command used is:\n{e.cmd}")
 
        sys.exit(1)
 
 
    return outFile


def queryExAC(variant:str)->str:

    """Query the ExAC database for one variant at a time.

    Parameters:
    variant       : str --- a string representation of a variant, in the form
                            chr-pos-ref-alt, such as, 1-931393-G-T.

    Returns:
    A string. a string representation of allele frequency exactly as in ExAC,
              such as '0.8439840567649248' or 'allele_freq_not_in_ExAC' when
              no such information in the database.
    """  

    ExAC_site="http://exac.hms.harvard.edu/rest/variant/variant/"

    url=ExAC_site+variant

    resp=requests.get(url)

    if resp.status_code !=200:

        print("Access to ExAC server failed for >variant<. Exit.")
        sys.exit(1)

    ExAC_info=resp.json()

    # if frequency not exists in database, let the user know
    if 'allele_freq' in ExAC_info:
        af=ExAC_info['allele_freq']
    else:
        af='allele_freq_not_in_ExAC'

    af_str=f"{af}" # turn into a string

    # a short delay to not choke hold the server
    # :)no need, they are very mean already :)
    time.sleep(0.01)


    return af_str


def parseInfo(info:str)->dict:

    """parse the 'INFO' field in VCF file into a dictionary."""

    info_dict={}

    temp=info.split(";")

    for x in temp:

       temp2=x.split("=")

       info_dict[temp2[0]]=temp2[1]

    return info_dict


def parseVepVcf(severity_dict, inFile:str, filter:bool=False)->list:

    """parse the VCF file from a vep run, one line at a time.

    Depth and read info for all the alleles will be parsed out from this file,
    along with the consequence and variant type information that were
    obtained by vep and stuck in INFO as a 'CSQ'.
 
    It will also query the ExAC database for allele frequecy.

    Parameters:
    severity_dict : dict  --- contains all the consequences and their
                              corresponding severity levels according to 
                              Ensembl. The lower the level, the more
                              severe.
    inFile        : str   --- the VCF file name
    filter        : bool  --- whether variants failing FILTER should be
                              ignored.
                              default: False

    Returns:
    A list     --- each item in the list is the information for one variant.
                   this list is used to output into a tsv file.
    """

    # a list holding all the info for output
    out=[]
    with open(inFile, 'r') as vcf:

        line=vcf.readline().rstrip("\n")

        # if not at the end of the file
        while line !='':

            # skip all the comment/info lines  
            if line.startswith("#"):

                line=vcf.readline().rstrip("\n")
                
                continue

            # split the fields for each variant in the VCF file
            field=line.split("\t")

            # ignore the ones that did not pass the filter, if filter==True
            if filter & (field[6] !='.') & (field[6].upper() !='PASS'):

                line=vcf.readline().rstrip("\n")

                continue
  
            # get all the alt alleles
            allele_alt_list=field[4].split(',')
 
            info_dict=parseInfo(field[7])

            # read depth at the site, all samples
            depth=info_dict['DP']
            
            # reference count
            count_ref=info_dict['RO']

            # get all the counts for the alt alleles
            count_alt=info_dict['AO']
            count_alt_list=count_alt.split(',')

            # match the alt count to the corresponding alt allele
            allele_count_alt_dict={}
            for i in range(len(allele_alt_list)):
                allele_count_alt_dict[allele_alt_list[i]]=count_alt_list[i]
            
            # calculate the percentages for all the alt alleles
            percentage_alt_list=[]
            for i in range(len(allele_alt_list)):
                
                percentage=float(count_alt_list[i])/float(depth)*100

                # keep all the decimals, to match to ExAC database
                percentage=f"{percentage}"

                percentage_alt_list+=[percentage]

            # turn into a string for output
            percentage_alt=','.join(percentage_alt_list)

            # this list holds all the info from ExAC for the variant
            ExAC_freq_list=[]

            # for each alt allele, query ExAC
            # this is silly. But bulk query always return a 405 error.
            # see limits in the accompanying README for more explanation.
            for alt in allele_alt_list:

                variant_t='-'.join([field[0], field[1], field[3], alt])                
                af_t=queryExAC('-'.join([field[0], field[1]
                                                       , field[3], alt]))

                ExAC_freq_list+=[af_t]

            # pick the largest frequence
            af=pickL(ExAC_freq_list)

            # select the most severe consequences from vep output
            CSQ_list=info_dict['CSQ'].split(',')

            symbol, allele, type, csq=selectMostSevere(severity_dict, CSQ_list)

            # calculate the percentage for ref
            percentage_ref=float(count_ref)/float(depth)*100
            percentage_ref=f"{percentage_ref}"

            # keep the first five fields in the VCF file in output
            origin_top="\t".join(field[0:5])

            line_t="\t".join([origin_top, symbol, depth, count_ref
                                        , percentage_ref
                                        , count_alt
                                        , percentage_alt
                                        , type, af, csq])

            out+=[line_t]

            # continue read one more line 
            line=vcf.readline().rstrip("\n")

    return out


def pickL(ExAC_freq_list:list)->str:

    """pick the largest float from the list, if any.
    otherwise, return 'allele_freq_not_in_ExAC'.
    """

    match=[]
    for f in ExAC_freq_list:

        if f != 'allele_freq_not_in_ExAC':
            match+=[float(f)] 

    if len(match)>0:
            
        af=max(match)
        af_str=f"{af}"

    else:
        af_str='allele_freq_not_in_ExAC'  


    return af_str


def selectMostSevere(severity_dict, CSQ0:list)->list:

    """Select the most severe consequence/effect from vep run.

    Parameters:
    severity_dict : dict  --- contains all the consequences and their
                              corresponding severity levels according to
                              Ensembl. The lower the level, the more
                              severe.

    CSQ0          : list  --- Each item in the list is the consequence string
                              from vep run per alt allele and transcript. It is
                              a '|' separated string, such as,
                              'CDK11B|SNV|intron_variant|T'.
 
                              There are four subitems there:

                              1st: gene symbol, if applicable, such as 'CDK11B'
                              2nd: variation_class, such as 'SNV'.
                              3rd: variant effect, such as 'intron_variant'
                              4th: alternate allele, such as 'C'

                              The 1st can be empty.

    Returns:  
    A list  -- The information about the most severe allele, in the order of
               gene symbol, alt allele, variant type and the most severe 
               consequence.
    """

    level=999 # a high number for no severity
    severity=''

    # loop through all the transcripts
    for trx in CSQ0:

        csq_trx_t=trx.split("|")

        symbol=csq_trx_t[0]
        type=csq_trx_t[1]
        allele=csq_trx_t[3]

        # loop through all the matches within a transcript
        if csq_trx_t[2]=='': # no consequence info
            continue

        # more than one consequence is possible per variant/transcript 
        for m in csq_trx_t[2].split('&'):

            if m in severity_dict:
                level_t=severity_dict[m]
            # in case there is a novel consequence, assign a high level
            # ** this does not bear any scientific significance **
            else:
                level_t=99

            if level>level_t:

                severity=m
                level=level_t

    return [symbol, allele, type, severity]


def setupSeverityTable(severityFile:str=None)->dict:

    """Build a dict to hold the consequence and its severity level.

       It will read in a two column tab delimited file. The first column is 
    the consequence and the second is the severity level (positive integer).
    The lower the level the more severe it is. The severity level is based
    on the information from Ensembl. The following is the web page:

         https://uswest.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences

       The default file is "myAnnotator.severity.txt". But one can supply
    their own via option '-s '.
 
       There can be comments in the file, starting with '#'. 

       *** The first line of the file must be a header. ***
       

    Parameters:
    severityFile      : str --- file name.

    Returns:
    A dict    --- the key is the consequence and the value is the level.
    """

    if severityFile==None:
        severityFile="myAnnotator.severity.txt"

    severity_level=[]
    severity_dict={}
    with open(severityFile, 'r') as severeFile:

        # skip the header
        line=severeFile.readline()

        line=severeFile.readline()

        while line !='':

            # skip the comment line 
            if line.startswith('#'):

                # next line
                line=severeFile.readline()

                continue
        
            temp=line.rstrip("\n").split("\t")

            severity_dict[temp[0]]=int(temp[1])
 
            line=severeFile.readline()
            # next line


    return severity_dict


def saveListToTabFile(anno:list, outFile:str):

    """Save the output into a tsv file.

    Parameters:
    anno          : list  --- Each item has all the output for a variant
    outFile       : str   --- output file name
    """

    # set up the file headers.
    header_list=['Chr', 'Pos', 'ID', 'Ref', 'Alt', 'Gene_symbol']
    header_list+=['Depth', 'Ref_Count', 'Ref_Percentage']
    header_list+=['Alt_Count', 'Alt_Percentage', 'Variant_Class']
    header_list+=['ExAC_freq']
    header_list+=['Effect']

    anno.insert(0, "\t".join(header_list))

    anno_str="\n".join(anno)
    with open(outFile, 'w') as saveList:

       saveList.write(anno_str)


if __name__=='__main__':

    """Program entry point"""
    main()

