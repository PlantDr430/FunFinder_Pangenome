#!/usr/bin/python3

'''
Analyses:
tmhmm
phobius
signalp
effectors
merops
dbcan
pfam
iprscan
antismash
'''

import os, re, sys, shutil, argparse, subprocess, itertools, random, textwrap, csv, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
from multiprocessing import Pool
from scipy.optimize import curve_fit
from collections import OrderedDict, defaultdict
from collections.abc import Iterable
from scipy.optimize import differential_evolution

currentdir = os.getcwd()
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -d directory -o output',
    description = '''    Data analysis of Funannotate and Orthofinder 
    outputs to Pangenome results.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-d',
    '--directory',
    required = True,
    help = 'Directory containing folders for processing',
    metavar=''
)
parser.add_argument(
    '-a',
    '--analyses',
    default = 'all',
    nargs='+',
    help = 'Analyses to compare [default: all]. Can put "none" or individual analyses separated '\
    'by commas',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required = True,
    help = 'New output folder. Do not use input directory as output directory',
    metavar=''
)
parser.add_argument(
    '-f_p',
    '--fischer_pvalue',
    default = 0.05,
    type = float,
    help = 'Pvalue cutoff for Fischer exact test (annotation bar charts) [default: 0.05]',
    metavar=''
)
parser.add_argument(
    '-k_p',
    '--kruskal_pvalue',
    default = 0.05,
    type = float,
    help = 'Pvalue cutoff for Kruskal-Wallis test (protein length boxplot) [default: 0.05]',
    metavar=''
)
parser.add_argument(
    '-b_p',
    '--benjamini_pvalue',
    default = 0.05,
    type = float,
    help = 'Pvalue cutoff for Benjamini-Hochberg in GO Enrichment Analysis [default: 0.05]',
    metavar=''
)
parser.add_argument(
    '-al',
    '--alpha',
    default = 0.05,
    type = float,
    help = 'Test-wise alpha for GO Enrichment Analysis [default: 0.05]',
    metavar=''
)
parser.add_argument(
    '-p',
    '--percent',
    default = 0.50,
    type = float,
    help = 'Percent of isolates within a cluster containing proteins with analysis hit, '\
    'to consider the cluster as such [default: 0.50]',
    metavar=''
)
parser.add_argument(
    '-c',
    '--cpus',
    type=int,
    default=1,
    help = 'Number of cores to use for multiprocessing of fluidity [default: 1]',
    metavar=''
)
parser.add_argument(
    '-max_sub',
    '--max_subsamples',
    type=int,
    default=50000,
    help = 'Max number of subsamples to run on N genomes sampled for fluidity. [default: 50000]',
    metavar=''
)
parser.add_argument(
    '--max_off',
    action='store_true',
    help = 'Turn off the max subsamples. This will cause the script sample ALL possible combinations'\
    'for N genomes',
)
parser.add_argument(
    '--NO_GO',
    action='store_true',
    help = 'Do not perform GO enrichment analysis if iprscan is in analyses [default: ON]',
)
parser.add_argument(
    '--EFFECTORP_PATH',
    help = 'Path to EffectorP.py if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--ENRICHMENT_PATH',
    help = 'Path to find_enrichment.py if not set in $PATH',
    metavar=''
)
parser.add_argument(
    '--GOBASIC_PATH',
    help = 'Path to go-basic.obo if not in current run directory',
    metavar=''
)
parser.add_argument(
    '--GOSLIM_PATH',
    help = 'Path to goslim_generic.obo if not in current run directory',
    metavar=''
)
args=parser.parse_args()

##### Preliminary set-up, make sure files, software, directories are present/created #####
def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

def walklevel(input_dir, level=1):
    input_dir = input_dir.rstrip(os.path.sep)
    if os.path.isdir(input_dir):
        assert os.path.isdir(input_dir)
        num_sep = input_dir.count(os.path.sep)
        for root, dirs, files in os.walk(input_dir):
            yield root, dirs, files
            num_sep_this = root.count(os.path.sep)
            if num_sep + level <= num_sep_this:
                del dirs[:]
    else:
        None

# Create output folder and define paths
input_dir = os.path.abspath(os.path.join(currentdir, args.directory))
input_directory = [os.path.abspath(x[0]) for x in walklevel(input_dir)]
input_directory.remove(os.path.abspath(input_dir))

if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out, 'species_results'))
    os.makedirs(os.path.join(args.out, 'working_directory', 'selection_results'))
result_dir = os.path.abspath(os.path.join(currentdir, args.out))
spec_dir = os.path.abspath(os.path.join(result_dir, 'species_results'))
work_dir = os.path.abspath(os.path.join(result_dir, 'working_directory'))
select_dir = os.path.abspath(os.path.join(work_dir, 'selection_results'))
if os.path.isdir(args.out):
    dirs = [os.path.join(args.out, 'species_results'), 
    os.path.join(args.out, 'working_directory', 'selection_results')]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)

for dir in os.listdir(input_dir):
    if os.path.isdir(os.path.join(input_dir, dir)):
        if not os.path.isdir(os.path.join(spec_dir, dir)):
            os.makedirs(os.path.join(spec_dir, dir))

# Determine input analyses
if args.analyses == 'all':
    analyses = ['tmhmm','phobius','signalp','effectors','merops','dbcan','pfam','iprscan','antismash']
elif args.analyses != 'all' and ''.join(args.analyses) != 'none':
    analyses = ''.join(args.analyses).lower().replace(' ','').split(',')
elif ''.join(args.analyses) == 'none':
    analyses = 'none'
else:
    print("ERROR: No analyses were provided, please provide analyses or 'none', -a or --analyses")
    sys.exit()

# Check to make sure files for specified analyses exist
if not analyses == 'none':
    iso_locus_prefix_dict = {}
    combined_fasta = os.path.abspath(os.path.join(work_dir, 'All.fasta'))
    if os.path.exists(combined_fasta):
        os.remove(combined_fasta) # clear old concatenated fasta if present, just in case
    for dir in input_directory:
        isolate = os.path.basename(dir)
        isolate_dir = os.path.abspath(os.path.join(spec_dir, isolate))
        if not os.path.isdir(os.path.join(isolate_dir, 'fasta_results')):
            os.makedirs(os.path.join(isolate_dir, 'fasta_results'))
        iso_fasta_dir = os.path.abspath(os.path.join(isolate_dir, 'fasta_results'))
        if not os.path.isdir(os.path.join(isolate_dir, 'work_dir')):
            os.makedirs(os.path.join(isolate_dir, 'work_dir'))
        iso_work_dir = os.path.abspath(os.path.join(isolate_dir, 'work_dir'))
        for a in analyses:
            found_file = [files.lower() for files in os.listdir(dir) if a in files.lower()] # look for file with analysis
            if a == 'antismash':
                if len(found_file) == 2: # There should be two antismash files that we will concatenate into one
                    anti_1 = os.path.abspath(os.path.join(dir, found_file[0]))
                    anti_2 = os.path.abspath(os.path.join(dir, found_file[1]))
                    antismash_file = os.path.abspath(os.path.join(iso_work_dir, 'antismash_annotations.txt'))
                    os.system('cat {} {} > {}'.format(anti_1,anti_2,antismash_file))
                else:
                    print('ERROR: Did not find both antismash.txt and antismashclusters.txt ',\
                    'file for isolate {} in {}'.format(isolate, dir))
                    sys.exit()
            elif a == 'effectors': # If specifying effectors, look for executable
                if args.EFFECTORP_PATH:
                    EFFECTORP = args.EFFECTORP_PATH
                else:
                    try:
                        if which_path('EffectorP.py'):
                            EFFECTORP = 'EffectorP.py'
                        elif which_path('EffectorP'):
                            EFFECTORP = 'EffectorP'
                        else:
                            raise
                    except:
                        print('ERROR: Effector analysis flag was given, but neither effector file nor EffectorP ',\
                        'executable found. Please provide pre-computed effector file or make sure parent directoy ,'\
                        'of EffectorP executable is located in $PATH')
                        sys.exit()
            elif a not in ''.join(found_file): # check to see if each analysis file is found for analyses in args
                print('ERROR: Did not find {} file for isolate {} in {}'.format(a, isolate, dir))
                sys.exit()
        fasta_ext = ['.fasta', '.fa', '.fas', '.fna']
        for ext in fasta_ext: # look for protein fasta file for each isolate
            found_file = [files for files in os.listdir(dir) if ext in files]
            if not found_file:
                print('ERROR: Did not find protein FASTA file for isolate {} in {}'.format(isolate, dir))
                sys.exit()
            with open(os.path.join(dir,''.join(found_file)), 'r') as prot_fasta,open(combined_fasta, 'a') as cat_fasta:
                contents = prot_fasta.readlines()
                cat_fasta.write(''.join(contents))
                iso_locus_prefix_dict[isolate] = (''.join(contents[0]).split('_',1)[0][1:])
            break

# Check to make sure Orthogroups file is present
for files in os.listdir(input_dir):
    if 'orthogroups' in files.lower():
        ortho_input = os.path.abspath(os.path.join(input_dir, files))
try:
    if ortho_input == None:
        raise
except:
    print('ERROR: Cannot find Orthogroups.txt file in {}'.format(input_dir))
    sys.exit()

# Look for scripts and files needed for GOEA
if 'iprscan' in analyses and not args.NO_GO:
    if args.ENRICHMENT_PATH:
        ENRICHMENT = args.ENRICHMENT_PATH
    else:
        try:
            if which_path('find_enrichment.py'):
                ENRICHMENT = 'find_enrichment.py'
            else:
                raise
        except:
            print('ERROR: find_enrichment.py not found, please make sure parent directory of ',\
            'find_enrichment.py is located in $PATH or provide path to executable in command ',\
            'with --ENRICHMENT_PATH')
            sys.exit()
    if args.GOBASIC_PATH:
        GOBASIC = args.GOBASIC_PATH
    else:
        try:
            basic=[files.lower() for files in os.listdir(currentdir) if 'go' in files.lower() and 'basic' in files.lower()]
            if not basic:
                raise
            else:
                GOBASIC = os.path.abspath(os.path.join(currentdir, ''.join(basic)))
        except:
            print('ERROR: Could not find go-basic.obo in {}, please make sure go-basic.obo is in ,'\
            'current run directory or provide path to file in command with --GOBASIC_PATH'.format(currentdir))
            sys.exit()
    if args.GOSLIM_PATH:
        GOSLIM = args.GOSLIM_PATH
    else:
        try:
            slim=[files.lower() for files in os.listdir(currentdir) if 'go' in files.lower() and 'slim' in files.lower()]
            if not slim:
                raise
            else:
                GOSLIM = os.path.abspath(os.path.join(currentdir, ''.join(slim)))
        except:
            print('ERROR: Could not find goslim-generic.obo in {}, please make sure goslim-generic.obo ,'\
            'is in current run directory or provide path to file in command with --GOSLIM_PATH'.format(currentdir))
            sys.exit()

########## Start main run ##########

def create_protein_fasta_dict(all_fasta_file):
    '''
    Create protein fasta dictionary of all proteins from combined fasta file.
    '''
    all_protein_dict = defaultdict(str)
    with open(combined_fasta, 'r') as all_proteins:
        name = ''
        for line in all_proteins:
            if line.startswith('>'):
                if ' ' in line:
                    name = line.split(' ')[0][1:]
                else:
                    name = line[1:-1]
                continue
            all_protein_dict[name]+=line.strip()
    return all_protein_dict

def create_ortho_dictionaries(orthogroups_text_file):
    '''
    Create multiple dictionaries from orthogroups text file.
    '''
    ortho_prot_dict = OrderedDict() # {Protein cluster ID : all proteins in cluster}
    ortho_iso_dict = OrderedDict() # {Protein cluster ID : all isolates that have proteins in cluster}
    with open(orthogroups_text_file, 'r') as infile:
        ortho_list = [item.strip() for item in sorted(infile)]
        for line in ortho_list:
            iso_list = []
            prot_list = []
            if ':' in line:
                cluster, genes = line.split(':')
            elif '\t' in line:
                cluster, genes = line.split('\t', 1)
            else:
                cluster, genes = line.split(' ', 1)
            for match in re.finditer(r'([^\s]+)', genes):
                isolate = match.group(0).split('_')[0]
                protein = match.group(0)
                # if not protein in all_protein_dict.keys():
                    # print(('Warning: Protein {} was found in ortho clusters but was not found in fasta file')
                    # .format(protein))
                iso_list.append(isolate)
                prot_list.append(protein)
            ortho_iso_dict[cluster] = list(set(iso_list))
            ortho_prot_dict[cluster] = prot_list

    # create pangenome dictionary for gene clusters
    count_dict = {} # temp dictionary for creation of pangenome categories
    pan_cluster_dict = {'Core': [], 'Accessory': [], 'Singleton': []} # {Pangenome category : list of clusters}
    iso_num = int(max([len(x) for x in ortho_iso_dict.values()]))
    for cluster in ortho_iso_dict.keys():
        count = len(ortho_iso_dict[cluster])
        if count == iso_num:
            pan_cluster_dict['Core'].append(cluster)
        elif iso_num > count > 1:
            pan_cluster_dict['Accessory'].append(cluster)
        elif count == 1:
            pan_cluster_dict['Singleton'].append(cluster)
        if count in count_dict.keys():
            count_dict[count] = count_dict[count] + 1
        else:
            count_dict[count] = 1
    sortedCountList = sorted(count_dict.items(), key=lambda i: -i[0]) # sort the count_dict to a list

    # create pangenome dictionary for proteins
    pan_prot_dict = {'Core': [], 'Accessory': [], 'Singleton': []} # {Pangenome category : list of proteins}
    for keys, values in pan_cluster_dict.items():
        if keys == 'Core':
            for x in values:
                pan_prot_dict['Core'].append(ortho_prot_dict[x])
        if keys == 'Accessory':
            for x in values:
                pan_prot_dict['Accessory'].append(ortho_prot_dict[x])
        if keys == 'Singleton':
            for x in values:
                pan_prot_dict['Singleton'].append(ortho_prot_dict[x])

    # Flatten out list of dictionary values (Faster this way then .get method)
    pan_prot_dict['Core'] = list(itertools.chain.from_iterable(pan_prot_dict['Core']))
    pan_prot_dict['Accessory'] = list(itertools.chain.from_iterable(pan_prot_dict['Accessory']))
    pan_prot_dict['Singleton'] = list(itertools.chain.from_iterable(pan_prot_dict['Singleton']))

    return ortho_prot_dict,ortho_iso_dict,pan_cluster_dict,pan_prot_dict,iso_num,sortedCountList

def create_pangenome_stats(sortCount,stats_text_out,stats_fig_out,iso_per_cluster,clusters_per_PAN,proteins_per_PAN):
    Core = 0 # Number of Core clusters in dataset
    Accessory = 0 # Number of Accessory clusters in dataset
    Singleton = 0 # Number of Singleton clusters in dataset
    Pangenome = 0 # Number of all clusters in dataset
    with open(stats_text_out, 'w') as stats_out:
        stats_out.write("Species Combinations\tNumber of gene clusters\n")
        for item in sortCount:
            stats_out.write('{}\t{}\n'.format(item[0], item[1]))
        Core = sortCount[0][1]
        Accessory = sum([x[1] for x in sortCount[1:-1]]) # all but first and last 
        Singleton = sortCount[-1][1]
        Pangenome = Core + Accessory + Singleton
        stats_out.write('Core\t{}\n'.format(str(Core)+'('+str(round((Core/Pangenome),4)*100)+'%)'))
        stats_out.write('Accessory\t{}\n'.format(str(Accessory)+'('+str(round((Accessory/Pangenome),4)*100)+'%)'))
        stats_out.write('Singleton\t{}\n'.format(str(Singleton)+'('+str(round((Singleton/Pangenome), 4)*100)+'%)'))
        stats_out.write('Total Pangenome\t{}\n'.format(Pangenome))

    # create figure
    df=pd.read_csv(pangenome_stats, delimiter='\t', header=0) # read is stats and create figure
    df.drop(df.tail(4).index,inplace=True) # drop bottom 4 rows of input stats file
    cluster_data=(np.array([int(x) for x in df['Number of gene clusters'].tolist()])) # get x-axis values
    fig, ax = plt.subplots(figsize=(10,4))
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', which='major', color='white')
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    c_bar = plt.bar(df['Species Combinations'], cluster_data, width=0.85, linewidth=1)
    for i in range(0, len(c_bar)): # loop through bars for specific colors
        if i == 0:
            c_bar[i].set_color('red') # core bar
            c_bar[i].set_edgecolor('black')
        elif 0 < i < len(c_bar)-1:
            c_bar[i].set_color('orange') # accessory bars
            c_bar[i].set_edgecolor('black')
        elif i == len(c_bar)-1:
            c_bar[i].set_color('yellow') # singleton bars
            c_bar[i].set_edgecolor('black')
    plt.xlim(-1,iso_num)
    plt.xlabel('Number of genomes')
    plt.ylabel('Number of protein clusters')
    red_patch = mpatches.Patch(color='red', label='Core')
    orange_patch = mpatches.Patch(color='orange', label='Accessory')
    yellow_patch = mpatches.Patch(color='yellow', label='Singleton')
    plt.legend(handles=[red_patch,orange_patch,yellow_patch], framealpha=1.0)
    plt.tight_layout()
    plt.savefig(stats_fig_out)
    plt.close()
    
    # print out some basic summary files 
    with open(iso_per_cluster, 'w') as iso_p_clus:
        for key, value in ortho_iso_dict.items():
            iso_p_clus.write(key+'\t'+' '.join(value)+'\t'+str(len(value))+'\n')
    with open(clusters_per_PAN, 'w') as clus_p_PAN:
        for key, value in pan_cluster_dict.items():
            clus_p_PAN.write(key+'\t'+' '.join(value)+'\n')
    with open(proteins_per_PAN, 'w') as prot_p_PAN:
        for key, value in pan_prot_dict.items():
            prot_p_PAN.write(key+'\t'+' '.join(value)+'\n')

    return Core, Accessory, Singleton, Pangenome

def exponential(x, a, b, c):
    return a * np.exp(b * x) + c

def neg_exponential(x, a, b, c):
    return a * np.exp(-b * x) + c

def sumOfSquaredError(parameterTuple, x_values, y_curve_values, func):
    warnings.filterwarnings("ignore") # do not print warnings by genetic algorithm
    val = func(x_values, *parameterTuple)
    return np.sum((y_curve_values - val) ** 2.0)

def generate_Initial_Parameters(x_values, y_curve_values, func):
    # min and max used for bounds
    maxX = max(x_values)
    minX = min(x_values)
    maxY = max(y_curve_values)
    minY = min(y_curve_values)
    maxXY = max(maxX, maxY)

    parameterBounds = []
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for a
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for b
    parameterBounds.append([-maxXY, maxXY]) # seach bounds for c
    # "seed" the numpy random number generator for repeatable results
    result = differential_evolution(sumOfSquaredError, parameterBounds, args=(x_values,y_curve_values, func), seed=3)
    return result.x

def subsample_multiprocess(combo_list):
    '''
    Takes portions of the full combo_list and runs them on separate threads for faster processing. 
    Calcualtes fluidity for each sample and returns list of fluidities.
    '''
    N = len(combo_list[0]) # get N from number of genomes present
    sample_process_list = []
    for sample in combo_list:
        pairs = tuple(itertools.combinations(sample,2))
        pair_fluidity_list = [pair_dict[tuple(sorted(p))] for p in pairs]
        sample_fluidity = (2/(N*(N-1)))*sum(pair_fluidity_list)
        sample_process_list.append(sample_fluidity)
    return sample_process_list

def flatten(lis):
    for item in lis:
     if isinstance(item, Iterable) and not isinstance(item, str):
         for x in flatten(item):
             yield x
     else:
         yield item

def create_pair_dictionary_for_fluidity():
    pair_dict = {} # {(Isolate1, Isolate2) : [ratio of sum(unique clusters)/sum(all clusters)]}
    for i in range(0, len(iso_list)):
        for x in range(0, len(iso_list)):
            if not iso_list[i] == iso_list[x]:
                pair = tuple(sorted([iso_list[i], iso_list[x]]))
                if not pair in pair_dict.keys():
                    cogs = {'Shared' : 0, 'Uk' : 0, 'Ul' : 0}
                    for v in ortho_iso_dict.values():
                        if pair[0] in v and pair[1] in v:
                            cogs['Shared'] += 1
                        elif pair[0] in v and pair[1] not in v:
                            cogs['Uk'] += 1
                        elif pair[0] not in v and pair[1] in v:
                            cogs['Ul'] += 1
                        else:
                            pass # don't need to count a cluster if both isolates are not present
                    unique_pair = cogs['Uk'] + cogs['Ul']
                    all_pair = (cogs['Uk'] + cogs['Shared']) + (cogs['Ul'] + cogs['Shared'])
                    pair_dict[pair] = unique_pair/all_pair
    return pair_dict

def determine_pangenome_fluidity(fluidity_figure, fluidity_results):
    '''
    Follows pangenome_fluidiy.py (https://github.com/PlantDr430/CSU_scripts/blob/master/pangenome_fluidity.py).
    '''
    # get fluidity and variance from total number of genomes (N)
    N = iso_num
    fluidity_list = [ratio for ratio in pair_dict.values()] # list of ratios 
    pan_fluidity = (2/(N*(N-1)))*sum(fluidity_list) # get fluidity from average of all ratios
    jack_samples = list(itertools.combinations(iso_list, N - 1)) # get list of all combos of N-1 from max num of genomes
    fluidity_i_list = []
    for sample in jack_samples:
        jack_pairs = tuple(itertools.combinations(sample,2)) # get all pairs from current jackknife sample
        jack_sample_fluidity = [pair_dict[tuple(sorted(p))] for p in jack_pairs] # get ratios from pair_dict
        fluidity_i = (2/((N-1)*(N-2)))*sum(jack_sample_fluidity) # calculate fluidity_i 
        fluidity_i_list.append(fluidity_i)
    fluidity_i_mean = np.mean(fluidity_i_list) # calculate fluidity_i_mean from all fluidity_i's
    pan_variance = ((N-1)/N)*sum([(i-fluidity_i_mean)**2 for i in fluidity_i_list]) # calculate variance

    # get fludiities to find varaiances for standard error of subsamples of randomly sampled 3 to N genomes
    sub_fluid_dict = {} # {N genomes sampled : [list of fluidities from subsamples]}
    for N in range(3, iso_num + 1):
        sub_fluid_dict[N] = []
        N_combos = list(itertools.combinations(iso_list, N))
        if args.max_off:
            combos = N_combos
        else:
            if len(N_combos) > args.max_subsamples:
                combos = random.choices(N_combos, k=args.max_subsamples)
                permutation_list.append(N) # create list of N genomes that hit max number of subsamples
            else:
                combos = N_combos
        if not len(N_combos) == 1:
            chunk = round(len(combos)/args.cpus)
            split_combos = [combos[i:i + chunk] for i in range(0, len(combos), chunk)]
            pool = Pool(processes=args.cpus)
            results = pool.imap(subsample_multiprocess, split_combos)
            sub_fluid_dict[N].append(results)
        else:
            last_run = subsample_multiprocess(N_combos)
            sub_fluid_dict[N].append(last_run)
        sub_fluid_dict[N]=list(flatten(sub_fluid_dict[N]))

    # create figure and print out results
    total_variance = []
    for i in range(3, iso_num + 1):
        if i in permutation_list:
            total_variance.append(np.var(sub_fluid_dict[i], ddof = 1) + pan_variance)
        else:
            total_variance.append(np.var(sub_fluid_dict[i]) + pan_variance)
    total_variance = np.array(total_variance)
    total_stderr = np.array([x**(1/2) for x in total_variance])
    y_fluidity_values = np.array([pan_fluidity for i in range(3, iso_num + 1)])
    x_labels = np.array([i for i in range(3, iso_num + 1)])
    stderr_bottom = np.array([(pan_fluidity - v) for v in total_stderr])
    stderr_top = np.array([(pan_fluidity + v) for v in total_stderr])
    fig, ax = plt.subplots()
    # try to determine best initial parameters for curve fittin
    try: # Still got problems sometimes with fitting curves, this temporary solution seems to be working
        geneticParameters_top = generate_Initial_Parameters(x_labels, stderr_top, exponential)
        geneticParameters_bottom = generate_Initial_Parameters(x_labels, stderr_bottom, exponential)
        popt_t, pcov = curve_fit(exponential, x_labels, stderr_top, geneticParameters_top, maxfev=10000)
        popt_b, pcov = curve_fit(exponential, x_labels, stderr_bottom, geneticParameters_bottom, maxfev=10000)
        if len(set(exponential(x_labels, *popt_t))) > 3 and len(set(exponential(x_labels, *popt_b))) > 3:
            plt.fill_between(x_labels, exponential(x_labels, *popt_t), exponential(x_labels, *popt_b), facecolor='blue', alpha=0.6)
            top_curve = exponential(x_labels, *popt_t)
            bottom_curve = exponential(x_labels, *popt_b)
        if len(set(exponential(x_labels, *popt_t))) <= 3:
            geneticParameters_top = generate_Initial_Parameters(x_labels, stderr_top, neg_exponential)
            popt_t, pcov = curve_fit(neg_exponential, x_labels, stderr_top, geneticParameters_top, maxfev=10000)
            plt.fill_between(x_labels, neg_exponential(x_labels, *popt_t), exponential(x_labels, *popt_b), facecolor='blue', alpha=0.6)
            top_curve = neg_exponential(x_labels, *popt_t)
            bottom_curve = exponential(x_labels, *popt_b)
        else:
            pass
        if len(set(exponential(x_labels, *popt_b))) <= 3:
            geneticParameters_bottom = generate_Initial_Parameters(x_labels, stderr_bottom, neg_exponential)
            popt_b, pcov = curve_fit(neg_exponential, x_labels, stderr_bottom, geneticParameters_bottom, maxfev=10000)
            plt.fill_between(x_labels, exponential(x_labels, *popt_t), neg_exponential(x_labels, *popt_b), facecolor='blue', alpha=0.6)
            top_curve = exponential(x_labels, *popt_t)
            bottom_curve = neg_exponential(x_labels, *popt_b)
        else:
            pass
    except:
        pass
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
    ax.xaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white', alpha=0.5)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    plt.plot(x_labels, y_fluidity_values, ls='--', lw=1, color='black') # plot y-values of fluidity
    plt.xticks(np.arange(x_labels[0], x_labels[len(x_labels)-1]+1, 1.0)) # make sure x interval is 1
    plt.xlim(x_labels[0], x_labels[len(x_labels)-1]) # adjust x limit so it starts with 3 at 0
    max_y = max(stderr_top)
    min_y = min(stderr_bottom)
    plt.ylim((min_y - min_y*0.15), (max_y + max_y*0.15))
    plt.xlabel('Number of genomes sampled')
    plt.ylabel('Fluidity, '+u'\u03C6')
    plt.tight_layout()
    plt.savefig(fluidity_figure)
    plt.close()

    with open(fluidity_results, 'w') as results: # print out fluidity results
        results.write('Genomes_Sampled\tFluidity\tTotal_Variance\tTotal_Stderr\tExponential_top\tExponential_bottom\n')
        r_out = []
        for i in range(0, iso_num-2):
            r_out.append([str(i+3), str(pan_fluidity), str(total_variance[i]), str(total_stderr[i]), 
            str(top_curve[i]), str(bottom_curve[i])])
        for line in r_out:
            results.write('\t'.join(line) + '\n')

def powerlaw(x, a, b, c):
    return a*(x**b) + c

def determine_pangenome_curve(curve_figure, curve_results):
    '''
    Determine pangenome curves. Core and pangenome counts of clusters for N randomly sampled genomes.
    '''
    isolate_dict = {} # {Isolate ID : list of protein clusters this isolate has}
    for id in iso_list:
        for key, value in ortho_iso_dict.items():
            if id in value:
                tmp_list = isolate_dict.get(id, [])
                isolate_dict[id] = tmp_list + [key]
        search_list = [id] # start with given isolate
        sub_list = iso_list.copy() # copy list of all isolates
        sub_list.remove(id) # remove start isolate from pool 
        random_sample = random.sample(sub_list, (int(iso_num) - 1)) # get random list from pool

        for i in range(0,(int(iso_num) - 2)):
            selection_out = os.path.abspath(os.path.join(select_dir, '{}-{}'.format(id,i+2)))
            if not os.path.exists(selection_out): # if slection files already exists, skip
                selection = random.choice(random_sample) # get a random isolate from random pool
                search_list.append(selection) # add choosen isolate to start isolate
                random_sample.remove(selection) # remove said choosen isolate from random pool to avoid duplicates

                sample_dict = {} # create pseudo orthogroups.txt for only selected isolates (from above)
                for key, value in ortho_iso_dict.items():
                    for x in search_list:
                        if x in value:
                            sample_list = sample_dict.get(key, [])
                            sample_dict[key] = sample_list + [x]

                with open(selection_out, 'w') as select_out: # print out selections to save memory
                    for key, value in sample_dict.items():
                        select_out.write(key + '\t' + '\t'.join(value) + '\n')

    curve_dict = {} # {IsolateID    N genomes sampled : [core genome count, pangenome count]}
    for files in os.listdir(select_dir): # loop through all selections and calculate pan/core genome values
        id, number = str(files).split('-')
        with open(os.path.join(select_dir, files)) as select:
            select_linearr = [item.strip() for item in sorted(select)]
            select_dict = {} # {Protein Cluster : List of isolates represented in cluster}
            for line in select_linearr:
                cluster, all_isolates = line.split('\t',1)
                for match in re.finditer(r'(\w+)', all_isolates):
                    isolates = match.group(0)
                    isolatesList = select_dict.get(cluster, [])
                    select_dict[cluster] = isolatesList + [isolates]
                if cluster not in select_dict.keys():
                    select_dict[cluster] = []
            sortedSelectList = sorted(select_dict.items(), key=lambda i: (-len(i[1]),i[0])) # sort list
            countDict = {} # {Number of Genomes : Number of clusters represented in this number of genomes}
            for key in select_dict.keys():
                count = len(select_dict[key])
                if count in countDict.keys():
                    countDict[count] = countDict[count] + 1
                else:
                    countDict[count] = 1
            sortedCount = sorted(countDict.items(), key=lambda i: -i[0])
            core = sortedCount[0][1]
            accessory = sum(x[1] for x in sortedCount[1:-1])
            singleton = sortedCount[-1][1]
            pangenome = core + accessory + singleton
            curve_dict[id + '\t' + number] = [str(core), str(pangenome)]

    # write out results from all N randomly sampled selections from above
    with open(curve_results, 'w') as curve_out:
        curve_out.write('Isolate\tGenomes\tCore\tPangenome\n')
        for key, value in curve_dict.items():
            curve_out.write(key + '\t' + '\t'.join(value) + '\n')
        for key, value in isolate_dict.items():
            curve_out.write(key + '\t' + '1' + '\t' + str(len(value)) + '\t' + str(len(value)) + '\n')
        curve_out.write('Total' + '\t' + str(iso_num) + '\t' + str(Core) + '\t' + str(Pangenome) + '\n')

    # create pangenome curve figure
    df=pd.read_csv(curve_results, delimiter='\t', header=0) # read curve results and create figure
    fig, ax = plt.subplots()
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', color='white', linestyle='--', alpha=0.3)
    ax.grid(which='major', linestyle='-', linewidth='1', color='white')
    ax.set_facecolor('gainsboro')
    x_data = df['Genomes']
    plt.scatter(x_data, df['Core'], c='red', label='Core')
    plt.scatter(x_data, df['Pangenome'], c='orange', label='Pangenome')
    x_power = np.array([float(x) for x in set(x_data.tolist())])
    core_avg = np.array(sorted(list(set(df.groupby(['Genomes']).Core.transform('mean').tolist())), reverse = True))
    pan_avg = np.array(sorted(list(set(df.groupby(['Genomes']).Pangenome.transform('mean').tolist()))))
    # guess intital parameters with a rough fit to a linear log fit. Tried using same method as fluidity but it failed
    rough_core = np.linalg.lstsq(np.stack((np.ones_like(x_power),np.log(x_power)),axis=1),np.log(core_avg),rcond=None)[0]
    p0_core = [np.exp(rough_core[0]), rough_core[1], 0]
    rough_pan = np.linalg.lstsq(np.stack((np.ones_like(x_power),np.log(x_power)),axis=1),np.log(pan_avg),rcond=None)[0]
    p0_pan = [np.exp(rough_pan[0]), rough_pan[1], 0]
    popt_core, pcov_core = curve_fit(powerlaw, x_power, core_avg, p0_core, maxfev=5000)
    popt_pan, pcov_pan = curve_fit(powerlaw, x_power, pan_avg, p0_pan, maxfev=5000)
    plt.plot(x_power, powerlaw(x_power, *popt_core), '-', color='black')
    plt.plot(x_power, powerlaw(x_power, *popt_pan), '-', color='black')
    plt.xlabel('Number of genomes sampled')
    plt.ylabel('Number of gene clusters')
    plt.legend(framealpha=1.0)
    plt.tight_layout()
    plt.savefig(curve_figure)
    plt.close()

def sig_letters(arg_value, p1v2, p1v3, p2v3):
    '''
    Creates list of letters to be placed in charts for significane differentiation. Currently only supports 
    graphs of at most 3 bars.
    '''
    test_labels = []
    test_labels.append('a')
    if p1v2 <= arg_value:
        test_labels.append('b')
    elif p1v2 > arg_value:
        test_labels.append('a')
    if p1v3 <= arg_value and p2v3 <= arg_value:
        test_labels.append('c')
    elif p1v3 <= arg_value and p2v3 > arg_value:
        test_labels.append('b')
    elif p1v3 > arg_value and p2v3 <= arg_value:
        test_labels.append('a')
    elif p1v3 > arg_value and p2v3 > arg_value:
        if p1v2 <= arg_value:
            test_labels.append('ab')
        elif p1v2 > arg_value:
            test_labels.append('a')
    return test_labels

def determine_pangenome_protein_lengths(length_figure, protein_lengths):
    '''
    Determine proteins lengths for all proteins per pangenome category and create boxplots.
    '''
    core_lengths = []
    accessory_lengths = []
    single_lengths = []
    for protein, sequence in all_protein_dict.items():
        if prot_pan_dict[protein] == 'Core':
            core_lengths.append(len(sequence))
        elif prot_pan_dict[protein] == 'Accessory':
            accessory_lengths.append(len(sequence))
        elif prot_pan_dict[protein] == 'Singleton':
            single_lengths.append(len(sequence))

    # write out tab-demilited file. Columns headers as pangenome categories, rows are protein lengths
    max_length = max(len(core_lengths),len(accessory_lengths),len(single_lengths))
    with open(protein_lengths, 'w') as prot_lengths:
        report = csv.writer(prot_lengths, delimiter = '\t', lineterminator='\n')
        report.writerow(['Core', 'Accessory', 'Singleton'])
        for i in range(0, max_length):
            row = []
            if i < len(core_lengths):
                row.append(core_lengths[i])
            else:
                row.append('')
            if i < len(accessory_lengths):
                row.append(accessory_lengths[i])
            else:
                row.append('')
            if i < len(single_lengths):
                row.append(single_lengths[i])
            else:
                row.append('')
            report.writerow(row)

    # create boxplots of protein lengths
    df=pd.read_csv(protein_lengths, delimiter='\t', header=0)
    core_list = df['Core'].values.tolist()
    acc_list = df['Accessory'].values.tolist()
    single_list = df['Singleton'].values.tolist()
    w1v2, p1v2 = stats.kruskal(core_list, acc_list, nan_policy='omit')
    w1v3, p1v3 = stats.kruskal(core_list, acc_list, nan_policy='omit')
    w2v3, p2v3 = stats.kruskal(core_list, acc_list, nan_policy='omit')
    kruskal_labels = sig_letters(args.kruskal_pvalue, p1v2, p1v3, p2v3)
    fig, ax = plt.subplots()
    lengths = np.column_stack((df['Core'], df['Accessory'], df['Singleton']))
    mask = ~np.isnan(lengths)
    filter_lengths = [d[m] for d, m in zip(lengths.T, mask.T)]
    ax.set_xticklabels(['Core', 'Accessory', 'Singleton'])
    ax.set_axisbelow(True)
    plt.minorticks_on()
    plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
    ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_facecolor('gainsboro')
    ax.axvline(x=1.5,color='white', linewidth='1') # insert line to separate 1st and 2nd boxplot
    ax.axvline(x=2.5,color='white', linewidth='1') # insert line to separate 2nd and 3rd boxplot
    bp = plt.boxplot(filter_lengths, 0, '', patch_artist=True) # don't plot outliers
    colors = ['red', 'orange', 'yellow']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    plt.setp(bp['boxes'], linewidth=2)
    plt.setp(bp['whiskers'], color='black', linewidth=2)
    plt.setp(bp['medians'], color='black', linewidth=2)
    plt.setp(bp['caps'], color='black', linewidth=2)
    pos = range(len(bp['boxes']))
    for tick, label in zip(pos, kruskal_labels): # loop through and place significance letters above median position
        median_pos = bp['medians'][tick].get_ydata()[0]
        label_pos = median_pos - (median_pos * 0.01)
        ax.text(pos[tick] + 1, label_pos, 
        label, ha='center', va='bottom', fontsize=12)
    plt.xlabel('Pangenome category')
    plt.ylabel('Protein length (aa)')
    plt.tight_layout()
    plt.savefig(length_figure)
    plt.close()

def parse_phobius(infile, secret_dict, trans_dict):
    with open(infile, 'r') as files:
        for line in files:
            if not 'PREDICTION' in line.upper():
                column = line.split('\t')
                if column[2] == 'Y': #column[1] == '0' and :
                    secret_dict[column[0]] = ['Secreted']
                elif int(column[1]) > 0: # and column[2] != 'Y':
                    trans_dict[column[0]] = [column[1]]
    return

def parse_tmhmm(infile, anno_dict):
    with open(infile, 'r') as files:
        for line in files:
            column = line.split('\t')
            tm_hit = re.search(r'(PredHel=)(\d+)', column[4]).group(2)
            if int(tm_hit) > 0:
                anno_dict[column[0]] = [tm_hit]
    return

def parse_signalp(infile, anno_dict):
    with open(infile, 'r') as files:
        for line in files:
            if not line.startswith('#'):
                column = line.split()
                if column[9] == 'Y' and column[11] != 'SignalP-TM':
                    anno_dict[column[0]] = ['Secreted']
    return

def parse_annotations(infile, anno_dict):
    with open(infile, 'r') as files:
        for line in files:
            column = line.split('\t')
            if 'go' in column[1].lower() or column[1].lower() == 'name':
                continue
            elif 'product' in column[1].lower():
                product_list = anno_dict.get(column[0], [])
                anno_dict[column[0]] = product_list + [column[2].strip()]
            else:
                anno_list = anno_dict.get(column[0], [])
                anno_dict[column[0]] = anno_list + [column[2].strip().split(':')[1]]
    return

def parse_antismash(infile, all_m, sm, back):
    with open(infile, 'r') as files:
        for line in files:
            column = line.split('\t')
            if column[1].lower() == 'product':
                back_list = back.get(column[0], [])
                back[column[0]] = back_list + [column[2].strip()]
            elif 'cluster' in column[2].lower():
                cluster_list = all_m.get(column[0], [])
                all_m[column[0]] = cluster_list + [column[2].strip().split(':')[1]]
            else:
                sm_list = sm.get(column[0], [])
                sm[column[0]] = sm_list + [column[2].strip().split(':')[0]]
    return

def parse_go(infile, anno_dict, anno_count_dict):
    go_count = {} # {GO-ID : number of times it occurs in current genome}
    filter_go_dict = {} # {Protein ID : list of GOs only present at least 5 times in current genome}
    with open(infile, 'r') as go_file:
        go_linearr = [item.strip() for item in go_file]
        for line in go_linearr: # loop through once to get counts and for all_go dictionary
            column = line.split('\t')
            if 'go' in column[1].lower():
                go = 'GO:'+column[2].split('|')[1]
                go_list = anno_dict.get(column[0], [])
                anno_dict[column[0]] = go_list + [go]
                anno_dict[column[0]] = list(set(anno_dict[column[0]])) # remove duplicates if any
                if go in go_count.keys():
                    go_count[go] = go_count[go] + 1
                else:
                    go_count[go] = 1
        # for GO enrichment analysis keep GO's that are only represented at least 5 times in the genome
        for line in go_linearr: # loop through a 2nd time for filter_go_dict (keep GO's >= 5)
            column = line.split('\t')
            if 'go' in column[1].lower():
                go = 'GO:'+column[2].split('|')[1]
                if go_count[go] >= 5:
                    filter_list = filter_go_dict.get(column[0], [])
                    filter_go_dict[column[0]] = filter_list + [go]
                    filter_go_dict[column[0]] = list(set(filter_go_dict[column[0]])) # remove duplicates if any
    for gene, GOs in filter_go_dict.items():
        for key, value in ortho_prot_dict.items():
            if gene in value:
                for g in range(0, len(GOs)):
                    dict_list = anno_count_dict.get(key+'_'+GOs[g], [])
                    anno_count_dict[key+'_'+GOs[g]] = dict_list + [isolate]
                    anno_count_dict[key+'_'+GOs[g]] = list(set(anno_count_dict[key+'_'+GOs[g]])) # remove duplicates
    return

def parse_effectors(infile, anno_dict):
    with open(infile, 'r') as files:
        for line in files:
            if '|' in line:
                column = line.strip().split('|')
                name = column[0]
                probability = column[1].split(':')[1]
                anno_dict[name] = [probability]
    return

def dict_to_fasta(outfile, dictionary):
    tmp_dict = {}
    for k, v in all_protein_dict.items():
        if k in dictionary.keys():
            tmp_dict[k] = v
    with open(outfile, 'w') as files:
        for k in tmp_dict.keys():
            files.write('>' + k + '\n')
            files.write(textwrap.fill(''.join(tmp_dict[k]),width=80) + '\n')
    return

def analyses_to_clusters(dictionary, dictionary_out, outfile1, outfile2):
    temp_dict = {}
    for gene in dictionary.keys():
        for key, value in ortho_prot_dict.items():
            if gene in value:
                gene_list = temp_dict.get(key, [])
                temp_dict[key] = gene_list + [gene]
                dict_list = dictionary_out.get(key, [])
                dictionary_out[key] = dict_list + [isolate]
                dictionary_out[key] = list(set(dictionary_out[key]))
    for key, value in temp_dict.items():
        outfile1.write(key+':\t' + '\t'.join(value) + '\n')
    core_cluster = set(temp_dict.keys()) & set(pan_cluster_dict['Core'])
    acc_cluster = set(temp_dict.keys()) & set(pan_cluster_dict['Accessory'])
    single_cluster = set(temp_dict.keys()) & set(pan_cluster_dict['Singleton'])
    core_prot = set(list(itertools.chain(*temp_dict.values()))) & set(pan_prot_dict['Core'])
    acc_prot = set(list(itertools.chain(*temp_dict.values()))) & set(pan_prot_dict['Accessory'])
    single_prot = set(list(itertools.chain(*temp_dict.values()))) & set(pan_prot_dict['Singleton'])
    outfile2.write('Core clusters:\t' + ' '.join(core_cluster) + '\n')
    outfile2.write('Number of core clusters:\t' + str(len(core_cluster))+ '\n')
    outfile2.write('Accessory clusters:\t' + ' '.join(acc_cluster) + '\n')
    outfile2.write('Number of accessory clusters:\t' + str(len(acc_cluster))+ '\n')
    outfile2.write('Singleton clusters:\t' + ' '.join(single_cluster) + '\n')
    outfile2.write('Number of singleton clusters:\t' + str(len(single_cluster))+ '\n')
    outfile2.write('Core proteins:\t' + ' '.join(core_prot) + '\n')
    outfile2.write('Number of core proteins:\t' + str(len(core_prot))+ '\n')
    outfile2.write('Accessory proteins:\t' + ' '.join(acc_prot) + '\n')
    outfile2.write('Number of accessory proteins:\t' + str(len(acc_prot)) + '\n')
    outfile2.write('Singleton proteins:\t' + ' '.join(single_prot) + '\n')
    outfile2.write('Number of accessory proteins:\t' + str(len(single_prot))+ '\n')
    return

def main_analyses_loop(loop_dir, isolate_dir, iso_fasta_dir, iso_work_dir):
    for cluster, proteins in complete_results.items(): # complete results = final excel file
        cluster_size = str(len(proteins))
        for k, v in pan_cluster_dict.items():
            if cluster in pan_cluster_dict['Core']:
                category = 'Core'
            elif cluster in pan_cluster_dict['Accessory']:
                category = 'Accessory'
            else:
                category = 'Singleton'
        for prot in proteins:
            id = prot[0].split('_')[0] # should be locus prefix tag
            prot.insert(0, id) # will end up as prot[3]
            prot.insert(0,category) # will end up as prot[2]
            prot.insert(0,cluster_size) # will end up as prot[1]
            prot.insert(0,cluster) # will end up as prot[0]
            # so far: complete results = {Cluster1 : [[cluster, cluster_size, PAN_category, locus_ID, protein_ID],[]],
            #                             Cluster2 : [[cluster, cluster_size, PAN_category, locus_ID, protein_ID],[]]}

    # annotation dictionaries
    gene_products = {} # {Protein ID : gene product name}
    phobius_s = {} # {Protein ID : 'Yes' for a phobius secreted hit}
    signalp = {} # {Protein ID : 'Yes' for a signalP secreted hit}
    phobius_t = {} # {Protein ID : number of transmembrane domains}
    tmhmm = {} # {Protein ID : number of transmembrane domains}
    effectors = {} # {Protein ID : effector probability score}
    merops = {} # {Protein ID : list of MEROPS IDs}
    dbcan = {} # {Protein ID : list of CAZy IDs}
    pfam = {} # {Protein ID : list of PFAM IDs}
    iprscan = {} # {Protein ID : list of IPR IDs}
    all_go = {} # {Protein ID : list GO terms}
    smcog = {} # {Protein ID : smCOG IDs}
    metabolites = {} # {Protein ID : antimash cluster ID protein is found in}
    backbone = {} # {Protein ID : name of metabolite backbone}
    for files in os.listdir(loop_dir): # loop through and parse annotation files
        if 'genes-products' in files.lower():
            product_file = os.path.abspath(os.path.join(loop_dir, files))
            parse_annotations(product_file, gene_products)
        if 'tmhmm' in analyses and 'tmhmm' in files.lower():
            tmhmm_file = os.path.abspath(os.path.join(loop_dir, files))
            parse_tmhmm(tmhmm_file, tmhmm)
        if 'phobius' in analyses and 'phobius' in files.lower():
            phobius_file = os.path.abspath(os.path.join(loop_dir, files))
            parse_phobius(phobius_file, phobius_s, phobius_t)
        if 'signalp' in analyses and 'signalp' in files.lower():
            signalp_file = os.path.abspath(os.path.join(loop_dir, files))
            parse_signalp(signalp_file, signalp)
        if 'merops' in analyses and 'merops' in files.lower():
            merops_file = os.path.abspath(os.path.join(loop_dir, files))
            parse_annotations(merops_file, merops)
        if 'dbcan' in analyses and 'dbcan' in files.lower():
            dbcan_file = os.path.abspath(os.path.join(loop_dir, files))
            parse_annotations(dbcan_file, dbcan)
        if 'iprscan' in analyses and 'iprscan' in files.lower():
            iprscan_file = os.path.abspath(os.path.join(loop_dir, files))
            parse_annotations(iprscan_file, iprscan)
            if not args.NO_GO:
                parse_go(iprscan_file, all_go, go_count_dict)
        if 'pfam' in analyses and 'pfam' in files.lower():
            pfam_file = os.path.abspath(os.path.join(loop_dir, files))
            parse_annotations(pfam_file, pfam)
    if 'antismash' in analyses: # Set out of for loop as we concatenated file earlier so no need to loop
        antismash_file = os.path.abspath(os.path.join(iso_work_dir, 'antismash_annotations.txt'))
        parse_antismash(antismash_file, metabolites, smcog, backbone)

    # Interactions of some dictionaries for combined dictionaries
    phobius_s_list = set(phobius_s.keys())
    phobius_t_list = set(phobius_t.keys())
    signalp_list = set(signalp.keys())
    tmhmm_list = set(tmhmm.keys())
    iprscan_list = set(iprscan.keys())
    pfam_list = set(pfam.keys())
    smcog_list = set(smcog.keys())
    metabolites_list = set(metabolites.keys())
    backbone_list = set(backbone.keys())
    secret_inter = phobius_s_list & signalp_list - tmhmm_list - phobius_t_list # add phobius_t as double check
    trans_inter = phobius_t_list & tmhmm_list
    conserve_inter = iprscan_list | pfam_list
    metabolite_inter = metabolites_list - smcog_list - backbone_list
    secretome = { i : ['Yes'] for i in secret_inter} # {Protein ID : 'Yes' for secretome}
    transmembrane = {i : ['Yes'] for i in trans_inter} # {Protein ID : 'Yes' for transmembrane}
    conserved = {i : ['Yes'] for i in conserve_inter} # {Protein ID : 'Yes' for conserved}
    # unknown metabolites are found in antismash clusters but don't have hits to smCOGs or backbones
    unknown_metabolites = {i : ['Yes'] for i in metabolite_inter} # {Protein ID : for unknown metabolite} 

    # Write out fasta files from dictionaries
    secretome_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'secretome.fasta'))
    dict_to_fasta(secretome_fasta, secretome)
    backbone_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'metabolite_backbones.fasta'))
    dict_to_fasta(backbone_fasta, backbone)
    smcog_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'metabolite_smcogs.fasta'))
    dict_to_fasta(smcog_fasta, smcog)
    unmetabolites_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'metabolite_unknowns.fasta'))
    dict_to_fasta(unmetabolites_fasta, unknown_metabolites)
    transmembrane_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'transmembrane.fasta'))
    dict_to_fasta(transmembrane_fasta, transmembrane)
    conserved_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'conserved_proteins.fasta'))
    dict_to_fasta(conserved_fasta, conserved)
    merops_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'merop_proteins.fasta'))
    dict_to_fasta(merops_fasta, merops)
    dbcan_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'cazy_proteins.fasta'))
    dict_to_fasta(dbcan_fasta, dbcan)

    # Run EffectorP from secretome file created above and create dictionary and output fasta
    if 'effectors' in analyses:
        effector_out = os.path.abspath(os.path.join(iso_work_dir, 'effector.out'))
        effector_log = os.path.abspath(os.path.join(iso_work_dir, 'effector.log'))
        if not os.path.exists(effector_out):
            with open(effector_log, 'w') as eff_log:
                subprocess.call([EFFECTORP, '-i', secretome_fasta, '-o', effector_out]
                , cwd=dir, stdout=eff_log, stderr=eff_log)
        parse_effectors(effector_out, effectors)
        effector_fasta = os.path.abspath(os.path.join(iso_fasta_dir, 'effectors.fasta'))
        dict_to_fasta(effector_fasta, effectors)

    for cluster, proteins in complete_results.items(): # finish creating complete results dictionary
        for prot in proteins:
            if iso_locus_prefix_dict[isolate] == prot[3]:
                if prot[4] in gene_products.keys():
                    prot.append(','.join(gene_products[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in pfam.keys():
                    prot.append(','.join(pfam[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in iprscan.keys():
                    prot.append(','.join(iprscan[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in all_go.keys():
                    prot.append(','.join(all_go[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in merops.keys():
                    prot.append(','.join(merops[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in dbcan.keys():
                    prot.append(','.join(dbcan[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in backbone.keys():
                    prot.append(','.join(backbone[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in smcog.keys():
                    prot.append(','.join(smcog[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in metabolites.keys():
                    prot.append(','.join(metabolites[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in signalp.keys():
                    prot.append(','.join(signalp[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in phobius_s.keys():
                    prot.append(','.join(phobius_s[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in phobius_t.keys():
                    prot.append(','.join(phobius_t[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in tmhmm.keys():
                    prot.append(','.join(tmhmm[prot[4]]))
                else:
                    prot.append(' ')
                if prot[4] in secretome.keys():
                    prot.append(','.join(secretome[prot[4]]))
                else:
                    prot.append('No')
                if prot[4] in transmembrane.keys():
                    prot.append(','.join(transmembrane[prot[4]]))
                else:
                    prot.append('No')
                if prot[4] in effectors.keys():
                    prot.append('Yes')
                    prot.append(''.join(effectors[prot[4]]))
                else:
                    prot.append('No')
                    prot.append('N/A')

    # print out results from individual analyses for each analysis and each isolate
    dictionary_list = ['metabolites','dbcan','transmembrane','secretome','effectors','conserved','merops']
    for dict in dictionary_list:
        dict_analysis_out = os.path.abspath(os.path.join(isolate_dir, dict+'_clusters_&_proteins.txt'))
        dict_overall_out = os.path.abspath(os.path.join(isolate_dir, 'overall_'+dict+'_results.txt'))
        with open(dict_analysis_out, 'w') as dict_out, open(dict_overall_out, 'w') as over_out:
            if dict == 'secretome':
                analyses_to_clusters(secretome, secretome_count_dict, dict_out, over_out)
            elif dict == 'effectors':
                analyses_to_clusters(effectors, effector_count_dict, dict_out, over_out)
            elif dict == 'transmembrane':
                analyses_to_clusters(transmembrane, transmembrane_count_dict, dict_out, over_out)
            elif dict == 'conserved':
                analyses_to_clusters(conserved, conserved_count_dict, dict_out, over_out)
            elif dict == 'metabolites':
                analyses_to_clusters(metabolites, metabolite_count_dict, dict_out, over_out)
            elif dict == 'merops':
                analyses_to_clusters(merops, merops_count_dict, dict_out, over_out)
            elif dict == 'dbcan':
                analyses_to_clusters(dbcan, dbcan_count_dict, dict_out, over_out)
    return

def perform_gene_ontology_enrichment():
    # Perform gene ontology enrichment (GOEA)
    Core_study = os.path.abspath(os.path.join(work_dir, 'core_study.txt'))
    Accessory_study = os.path.abspath(os.path.join(work_dir, 'accessory_study.txt'))
    Singleton_study = os.path.abspath(os.path.join(work_dir, 'singleton_study.txt'))
    population = os.path.abspath(os.path.join(work_dir, 'population.txt'))
    go_assoc = os.path.abspath(os.path.join(work_dir, 'go_associations.txt'))
    total_core_go = []
    total_acc_go = []
    total_single_go = []
    association_dict = {} # {Protein cluster : list of GO terms in cluster}
    for cluster, isolates in go_count_dict.items():
        clus = cluster.split('_')[0]
        go = cluster.split('_')[1].strip()
        # if GO term contains a number of isolates >= args.percent for isolates in a cluster, that cluster has GO term
        if len(isolates)/len(ortho_iso_dict[clus]) >= args.percent:
            association_list = association_dict.get(clus, [])
            association_dict[clus] = association_list + [go]
            association_dict[clus] = list(set(association_dict[clus]))
            if cluster_pan_dict[clus] == 'Core':
                total_core_go.append(clus)
            elif cluster_pan_dict[clus] == 'Accessory':
                total_acc_go.append(clus)
            elif cluster_pan_dict[clus] == 'Singleton':
                total_single_go.append(clus)
    total_core_go = sorted(list(set(total_core_go)))
    total_acc_go = sorted(list(set(total_acc_go)))
    total_single_go = sorted(list(set(total_single_go)))
    with open(Core_study, 'w') as c_study:
        c_study.write('\n'.join(total_core_go))
    with open(Accessory_study, 'w') as a_study:
        a_study.write('\n'.join(total_acc_go))
    with open(Singleton_study, 'w') as s_study:
        s_study.write('\n'.join(total_single_go))
    with open(population, 'w') as pop_out:
        for k in ortho_iso_dict.keys():
            pop_out.write(k + '\n')
    with open(go_assoc, 'w') as assoc_out:
        for k, v in sorted(association_dict.items()):
            assoc_out.write(k + '\t' + ';'.join(v) + '\n')

    study_list = [Core_study, Accessory_study, Singleton_study]
    for study in study_list:
        category = os.path.basename(study).split('_')[0]
        goea_log = os.path.abspath(os.path.join(work_dir, 'GOEA_{}.log'.format(category)))
        goea_out = os.path.abspath(os.path.join(work_dir, 'GOEA_{}.tsv'.format(category)))
        with open(goea_log, 'w') as logfile:#, with open(study, 'r') as study_file:
            subprocess.call([ENRICHMENT, study, population, go_assoc, '--method=fdr_bh', 
            '--pval_field=fdr_bh', '--alpha', str(args.alpha), '--pval', str(args.benjamini_pvalue), 
            '--obo', GOBASIC, '--goslim', GOSLIM, '--no_propagate_counts', '--outfile', goea_out
            ], cwd=currentdir, stdout=logfile, stderr=logfile)
    
    category_list = ['Core', 'Accessory', 'Singleton']
    type_dict = {'BP':'Biological Process', 'MF': 'Molecular Function', 'CC':'Cellular Component'}
    for files in os.listdir(work_dir):
        for cat in category_list:
            if cat.lower()+'.tsv' in files:
                fpath = os.path.abspath(os.path.join(work_dir, files))
                df = pd.read_csv(fpath, delimiter='\t', header=0)
                df.drop(df[df['enrichment'] == 'p'].index, inplace=True) # get only rows that are enriched
                for type in type_dict.keys():
                    goea_figure = os.path.abspath(os.path.join(result_dir, 'GOEA_{}_{}.png'.format(cat,type)))
                    GO_type = df['NS'] == type
                    type_df = df[GO_type] # created new dataframes based on current type
                    if not type_df.empty:
                        log_Pvalues = np.flip(np.array([[(-np.log10(float(x)))] for x in type_df['p_fdr_bh'].tolist()]))
                        len_P = len(log_Pvalues)
                        go_ticks = [g for g in type_df['# GO'].tolist()] # get GO names
                        all_y_ticks = [n for n in type_df['name'].tolist()] # get function names
                        y_ticks = []
                        for tick in all_y_ticks:
                            if ',' in tick:
                                short_tick = tick.split(',')[0] # only use first part of name if commas are present
                                y_ticks.append(short_tick)
                            else:
                                y_ticks.append(tick)
                        # get best guesses for figure size and aspect ratio based on number of y-values
                        fig_size = (8+(len(y_ticks)*0.021), len(y_ticks)*(.4+(1/(len(y_ticks)+1))))
                        aspect_ratio = len(y_ticks)*(0.3+(1/(len(y_ticks)+10)))
                        bar_shrink = 1/(2+((len(y_ticks)/100)*len(y_ticks)*0.03))
                        bar_aspect =  2.5+(len(y_ticks)/5)

                        fig, ax= plt.subplots(figsize=fig_size)
                        im = plt.imshow(log_Pvalues, cmap='Blues', aspect=aspect_ratio, extent=(0,len_P,0,len_P))
                        ax.get_xaxis().set_visible(False)
                        ax.set_yticks(ticks=np.arange(0.5, len(y_ticks)))
                        ax.set_yticklabels(y_ticks, fontsize=12)
                        x_value = ax.get_xlim()
                        for i in range(len(y_ticks)):
                            text = ax.text(np.mean(x_value), i+0.5, go_ticks[i], ha='center', va='center', fontsize=10)
                        if len(y_ticks) <= 4: # slightly different colorbar for plots with low number of y-values
                            cbar = plt.colorbar(im, pad=0.15, shrink=bar_shrink,aspect=bar_aspect)
                        else:
                            # scalar_ticks = np.linspace(0,np.round(max(log_Pvalues)), 4, dtype='int')
                            # cbar = plt.colorbar(shrink=bar_shrink, pad=0.15, aspect=bar_aspect, ticks=scalar_ticks)
                            cbar = plt.colorbar(shrink=bar_shrink, pad=0.15, aspect=bar_aspect)
                        cbar.ax.tick_params(labelsize=12)
                        cbar.ax.invert_yaxis()
                        cbar.ax.set_title('-log10(P-value)', fontsize=12)
                        plt.title('{} {}'.format(cat,type_dict[type]), fontsize=16, pad=10)
                        plt.tight_layout()
                        plt.savefig(goea_figure)
                        plt.close()

def analyses_report(input_file, count_dictionary):
    core_list = []
    acc_list = []
    single_list = []
    for i, line in enumerate(input_file):
        if i == 0:
            category, core_clusters = line.split(':\t')
            for core_hit in re.finditer(r'(\w+\d+)', core_clusters):
                if core_hit.group(0) in count_dictionary.keys():
                    if (len(count_dictionary[core_hit.group(0)])/
                    len(ortho_iso_dict[core_hit.group(0)])) >= args.percent:
                        core_list.append(core_hit.group(0))
        elif i == 2:
            category, acc_clusters = line.split(':\t')
            for acc_hit in re.finditer(r'(\w+\d+)', acc_clusters):
                if acc_hit.group(0) in count_dictionary.keys():
                    if (len(count_dictionary[acc_hit.group(0)])/
                    len(ortho_iso_dict[acc_hit.group(0)])) >= args.percent:
                        acc_list.append(acc_hit.group(0))
        elif i == 4:
            category, single_clusters = line.split(':\t')
            for single_hit in re.finditer(r'(\w+\d+)', single_clusters):
                if single_hit.group(0) in count_dictionary.keys():
                    # if (len(count_dictionary[single_hit.group(0)])/
                    # len(ortho_iso_dict[single_hit.group(0)])) >= args.percent:
                    single_list.append(single_hit.group(0))
    return core_list, acc_list, single_list

def write_final_results_and_figure(loop_dir):
    dictionary_list = ['metabolites','dbcan','transmembrane','secretome','effectors','conserved','merops']
    for dict in dictionary_list:
        isolate = os.path.basename(loop_dir)
        isolate_dir = os.path.abspath(os.path.join(spec_dir, isolate))
        overall_results = os.path.abspath(os.path.join(isolate_dir, 'overall_'+dict+'_results.txt'))

        # read in results from main analyses loop
        with open(overall_results, 'r') as overall_analysis:
            if dict == 'secretome':
                secretome_results = analyses_report(overall_analysis, secretome_count_dict)
                total_core_secretome = secretome_results[0]
                total_acc_secretome = secretome_results[1]
                total_single_secretome = secretome_results[2]
            elif dict == 'effectors':
                effector_results = analyses_report(overall_analysis, effector_count_dict)
                total_core_effectors = effector_results[0]
                total_acc_effectors = effector_results[1]
                total_single_effectors = effector_results[2]
            elif dict == 'transmembrane':
                transmembrane_results = analyses_report(overall_analysis, transmembrane_count_dict)
                total_core_transmembrane = transmembrane_results[0]
                total_acc_transmembrane = transmembrane_results[1]
                total_single_transmembrane = transmembrane_results[2]
            elif dict == 'conserved':
                conserved_results = analyses_report(overall_analysis, conserved_count_dict)
                total_core_conserved = conserved_results[0]
                total_acc_conserved = conserved_results[1]
                total_single_conserved = conserved_results[2]
            elif dict == 'metabolites':
                metabolite_results = analyses_report(overall_analysis, metabolite_count_dict)
                total_core_metabolites = metabolite_results[0]
                total_acc_metabolites = metabolite_results[1]
                total_single_metabolites = metabolite_results[2]
            elif dict == 'merops':
                merops_results = analyses_report(overall_analysis, merops_count_dict)
                total_core_merops = merops_results[0]
                total_acc_merops = merops_results[1]
                total_single_merops = merops_results[2]
            elif dict == 'dbcan':
                dbcan_results = analyses_report(overall_analysis, dbcan_count_dict)
                total_core_dbcan = dbcan_results[0]
                total_acc_dbcan = dbcan_results[1]
                total_single_dbcan = dbcan_results[2]

    # write out final results
    final_results = os.path.abspath(os.path.join(result_dir, 'Pangenome_analyses.txt'))
    with open(final_results, 'w') as results:
        results.write('Total core transmembrane clusters:\t' + str(len(set(total_core_transmembrane))) + '\n')
        results.write('Total accessory transmembrane clusters:\t' + str(len(set(total_acc_transmembrane))) + '\n')
        results.write('Total singleton transmembrane clusters:\t' + str(len(set(total_single_transmembrane))) + '\n')
        results.write('Total core secretome clusters:\t' + str(len(set(total_core_secretome))) + '\n')
        results.write('Total accessory secretome clusters:\t' + str(len(set(total_acc_secretome))) + '\n')
        results.write('Total singleton secretome clusters:\t' + str(len(set(total_single_secretome))) + '\n')
        results.write('Total core effector clusters:\t' + str(len(set(total_core_effectors))) + '\n')
        results.write('Total accessory effector clusters:\t' + str(len(set(total_acc_effectors))) + '\n')
        results.write('Total singleton effector clusters:\t' + str(len(set(total_single_effectors))) + '\n')
        results.write('Total core 2nd metabolite clusters:\t' + str(len(set(total_core_metabolites))) + '\n')
        results.write('Total accessory 2nd metabolite clusters:\t' + str(len(set(total_acc_metabolites))) + '\n')
        results.write('Total singleton 2nd metabolite clusters:\t' + str(len(set(total_single_metabolites))) + '\n')
        results.write('Total core clusters with conserved protein domains:\t' + str(len(set(total_core_conserved))) + '\n')
        results.write('Total accessory clusters with conserved protein domains:\t' + str(len(set(total_acc_conserved))) + '\n')
        results.write('Total singleton clusters with conserved protein domains:\t' + str(len(set(total_single_conserved))) + '\n')
        results.write('Total core MEROP proteins:\t' + str(len(set(total_core_merops))) + '\n')
        results.write('Total accessory MEROP proteins:\t' + str(len(set(total_acc_merops))) + '\n')
        results.write('Total singleton MEROP proteins:\t' + str(len(set(total_single_merops))) + '\n')
        results.write('Total core CAZy proteins:\t' + str(len(set(total_core_dbcan))) + '\n')
        results.write('Total accessory CAZy proteins:\t' + str(len(set(total_acc_dbcan))) + '\n')
        results.write('Total singleton CAZy proteins:\t' + str(len(set(total_single_dbcan))) + '\n')
        results.write('Total core clusters:\t' + str(len(set(pan_cluster_dict['Core']))) + '\n')
        results.write('Total accessory clusters:\t' + str(len(set(pan_cluster_dict['Accessory']))) + '\n')
        results.write('Total singleton clusters:\t' + str(len(set(pan_cluster_dict['Singleton']))) + '\n')

    # Create bar graphs for all analyses from final results
    final_analyses = ['metabolite','CAZy','transmembrane','secretome','effector','conserved','MEROP']
    for a in final_analyses:
        analysis_graph = os.path.abspath(os.path.join(result_dir, a+'_pangenome_bar.png'))
        analysis_data = [] # percents of analysis out of PANcategory for bar graphs
        f_list = [] # list of values for fischer exact tests
        with open(final_results, 'r') as data_file:
            for line in data_file:
                header, prot_count = line.split('\t')
                if a in header and 'core' in header:
                    analysis_data.append((int(prot_count.strip())/Core)*100)
                    f_list.append([int(prot_count.strip())])
                elif a in header and 'accessory' in header:
                    analysis_data.append((int(prot_count.strip())/Accessory)*100)
                    f_list[0].append(int(prot_count.strip()))
                elif a in header and 'singleton' in header:
                    analysis_data.append((int(prot_count.strip())/Singleton)*100)
                    f_list[0].append(int(prot_count.strip()))
                    f_list.append([Core, Accessory, Singleton])
        one_v_two = [[f_list[0][0], f_list[0][1]],[f_list[1][0],f_list[1][1]]]
        one_v_three = [[f_list[0][0], f_list[0][2]],[f_list[1][0],f_list[1][2]]]
        two_v_three = [[f_list[0][1], f_list[0][2]],[f_list[1][1],f_list[1][2]]]
        o1v2, p1v2 = stats.fisher_exact(one_v_two)
        o1v3, p1v3 = stats.fisher_exact(one_v_three)
        o2v3, p2v3 = stats.fisher_exact(two_v_three)
        fischer_labels = sig_letters(args.fischer_pvalue, p1v2, p1v3, p2v3)
        label = ['Core','Accessory','Singleton']
        fig, ax = plt.subplots()
        ax.set_axisbelow(True)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())
        plt.minorticks_on()
        plt.grid(which='minor', axis='y', color='white', linestyle='--', alpha=0.3)
        ax.yaxis.grid(True, linestyle='-', linewidth='1', which='major', color='white')
        ax.tick_params(axis='x', which='minor', bottom=False)
        ax.set_xticklabels(['Core', 'Accessory', 'Singleton'])
        ax.axvline(x=0.5,color='white', linewidth='1')
        ax.axvline(x=1.5,color='white', linewidth='1')
        ax.set_facecolor('gainsboro')
        plt.bar(label, analysis_data, width=0.85, linewidth=2, color=['red','orange','yellow'], edgecolor='black')
        rects = ax.patches
        (y_bottom, y_top) = ax.get_ylim()
        y_height = y_top - y_bottom
        for rect, label in zip(rects, fischer_labels):
            height = rect.get_height()
            label_position = height + (y_height * 0.01)
            ax.text(rect.get_x() + rect.get_width() /2, label_position, 
            label, ha='center', va='bottom', fontsize = 12)
        plt.ylim(0,y_height*1.03)
        plt.xlabel('Pangenome categories')
        if a == 'metabolite':
            plt.ylabel(u'2\u00B0 metabolite protein clusters (%)')
        elif a == 'transmembrane':
            plt.ylabel('Transmembrane protein clusters (%)')
        elif a == 'secretome':
            plt.ylabel('Secreted protein clusters (%)')
        elif a == 'effector':
            plt.ylabel('Predicted effector protein clusters (%)')
        elif a == 'conserved':
            plt.ylabel('Protein clusters with conserved protein domains (%)')
        elif a == 'MEROP':
            plt.ylabel('MEROP protein clusters (%)')
        elif a == 'CAZy':
            plt.ylabel('CAZy protein clusters (%)')
        plt.tight_layout()
        plt.savefig(analysis_graph)
        plt.close()

def create_final_excel_file(excel_report):
    '''
    Combine results of analyses for all isolates into final pangenome report.
    '''
    with open(excel_report, 'w') as excel_out:
        excel = csv.writer(excel_out, delimiter = '\t', lineterminator='\n')
        excel.writerow(['Pangenome_cluster', 'Pangenome_cluster_size', 'Pangenome_category', 'Protein_id',
        'Gene_product', 'PFAM_domains', 'Iprscan_domains', 'Gene_ontology', 'MEROPS', 'CAZy', '2nd_metabolite_backbone', 
        'smCOGs', 'Antismash_cluster', 'SignalP', 'Phobius_secreted', 'Phobius_transmembrane_domains', 'TMHMM', 
        'Secretome', 'Transmembrane', 'EffecorP_candidate', 'EffectorP_score'])
        for cluster, proteins in complete_results.items():
            for prot in proteins:
                remove_locus = prot.pop(3) # remove locus_prefix as it was just used for indexing previously
                excel.writerow(prot)
    return


if __name__ == "__main__":
    #### dictionary creation ####
    print('Creating pangenome dictionaries and lists from orthogroups')
    all_protein_dict = create_protein_fasta_dict(combined_fasta)
    ortho_dictionaries = create_ortho_dictionaries(ortho_input)
    ortho_prot_dict = ortho_dictionaries[0] # {Protein cluster ID : all proteins in cluster}
    ortho_iso_dict = ortho_dictionaries[1] # {Protein cluster ID : all isolates that have proteins in cluster}
    pan_cluster_dict = ortho_dictionaries[2] # {Pangenome category : list of clusters}
    pan_prot_dict = ortho_dictionaries[3] # {Pangenome category : list of proteins}
    iso_num = ortho_dictionaries[4] # number of genomes in sample
    sortCount = ortho_dictionaries[5] # sortedCountList for pangenome stats
    iso_list = list(set(itertools.chain.from_iterable([v for v in ortho_iso_dict.values() if len(v) == iso_num])))

    # flip pan_cluster_dict for alternative searching {Cluster : Pangenome category}
    cluster_pan_dict = {}
    for pan, cluster in pan_cluster_dict.items():
        for clus in cluster:
            cluster_pan_dict[clus] = pan
    # flip pan_prot_dict for alternative searching {Protein : Cluster it is in}
    prot_pan_dict = {}
    for k, v in pan_prot_dict.items():
        for prot in v:
            prot_pan_dict[prot] = k

    print('Working on general pangenome reports from orthogroups')

    #### pangenome stats ####
    print('Writing general pangenome stats')
    pangenome_stats = os.path.abspath(os.path.join(result_dir, 'Pangenome_stats.txt'))
    stats_figure = os.path.abspath(os.path.join(result_dir, 'Pangenome_stats.png'))
    isolates_per_cluster = os.path.abspath(os.path.join(result_dir, 'Isolates_per_cluster.txt'))
    cluters_per_PANcategory = os.path.abspath(os.path.join(result_dir, 'Cluters_per_PANcategory.txt'))
    proteins_per_PANcategory = os.path.abspath(os.path.join(result_dir, 'Proteins_per_PANcategory.txt'))
    pan_stats_results = create_pangenome_stats(sortCount, pangenome_stats, stats_figure, isolates_per_cluster, 
                        cluters_per_PANcategory, proteins_per_PANcategory)
    Core = pan_stats_results[0]
    Accessory = pan_stats_results[1]
    Singleton = pan_stats_results[2]
    Pangenome = pan_stats_results[3]

    #### pangenome fluidity ####
    print('Determning pangenome fluidity')
    pangenome_fluidity_results = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.txt'))
    pangenome_fluidity_fig = os.path.abspath(os.path.join(result_dir, 'Pangenome_fluidity.png'))
    permutation_list = []
    pair_dict = create_pair_dictionary_for_fluidity()
    determine_pangenome_fluidity(pangenome_fluidity_fig, pangenome_fluidity_results)

    ##### pangenome curve ####
    print('Determning pangenome curve')
    pangenome_curve_results = os.path.abspath(os.path.join(result_dir, 'Pangenome_curve.txt'))
    pangenome_curve_figure = os.path.abspath(os.path.join(result_dir, 'Pangenome_curve.png'))
    determine_pangenome_curve(pangenome_curve_figure, pangenome_curve_results)

    #### pangenome protein lengths ####
    print('Determning proteins lengths for all proteins based on pangenome category')
    protein_lengths_results = os.path.abspath(os.path.join(result_dir, 'Proteins_lengths.txt'))
    protein_length_figure = os.path.abspath(os.path.join(result_dir, 'Protein_lengths.png'))
    determine_pangenome_protein_lengths(protein_length_figure, protein_lengths_results)

    # If analyses equal 'none' finish run, else start working through each species
    if analyses == 'none':
        print('Finished basic pangenome results')
    else:
        print('Finished basic pangenome results. Working on combining analyses for each isolate',\
        'with pangenome information for {} isolates.'.format(iso_num))

    #### loop for analyses ####
    secretome_count_dict = {} # {Protein cluster : list of isolates classifying cluster as secretome}
    effector_count_dict = {} # {Protein cluster : list of isolates classifying cluster as effector}
    metabolite_count_dict = {} # {Protein cluster : list of isolates classifying cluster as metabolite}
    transmembrane_count_dict = {} # {Protein cluster : list of isolates classifying cluster as transmembrane}
    conserved_count_dict = {} # {Protein cluster : list of isolates classifying cluster as conserved}
    merops_count_dict = {} # {Protein cluster : list of isolates classifying cluster as merops}
    dbcan_count_dict = {} # {Protein cluster : list of isolates classifying cluster as dbcan}
    go_count_dict = {} # {Protein cluster_GO-ID : list of isolates containing this GO-ID in this cluster}
    # complete_results = {Protein cluster : nested list of all proteins in cluster}
    complete_results = {cluster : [[prot] for prot in ortho_prot_dict[cluster]] for cluster in ortho_prot_dict.keys()}
    isolate_count = 0 
    for dir in input_directory:
        isolate_count += 1
        isolate = os.path.basename(dir)
        isolate_dir = os.path.abspath(os.path.join(spec_dir, isolate))
        iso_fasta_dir = os.path.abspath(os.path.join(isolate_dir, 'fasta_results'))
        iso_work_dir = os.path.abspath(os.path.join(isolate_dir, 'work_dir'))
        print('Working on {} {}/{} isolates'.format(isolate, isolate_count, iso_num))
        main_analyses_loop(dir, isolate_dir, iso_fasta_dir, iso_work_dir)

    #### Perform gene ontology enrichment ####
    if 'iprscan' in analyses and not args.NO_GO:
        print('Performing Gene Ontology Enrichment Analysis (GOEA)')
        perform_gene_ontology_enrichment()

    #### write final analyses results from finished loop results####
    print('Writing out final results from analyses and creating bar graphs')
    for dir in input_directory:
        write_final_results_and_figure(dir)

    #### Create final excel report ####
    print('Creating final excel spreadsheet of pangenomic information')
    pangenome_excel_results = os.path.abspath(os.path.join(result_dir, 'All_pangenome_info.txt'))
    create_final_excel_file(pangenome_excel_results)

