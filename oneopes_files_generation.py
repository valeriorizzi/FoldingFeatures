import argparse
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import os
import sys
import datetime
import shutil

error='--- ERROR: %s \n'

parser = argparse.ArgumentParser(
    prog = "oneopes_files_generation",
    description = "Generates OneOPES-ready files from plumed_final_noduplicate.dat, COLVAR_diff and COLVAR_solvation.",
    epilog = "Thanks for using %(prog)s! ",
)

def read_colvar(name: str):
    """
    Read the colvar files of the folded and unfolded directories
    Returns a dataframe with some statistics on loaded CVs (except time)
    """
    # Reads COLVAR file
    filenamesF = args_path_folded + "/" + name
    filenamesU = args_path_unfolded + "/" + name

    header = pd.read_csv(filenamesF, sep=r'\s+', nrows=0).columns.tolist()
    values_to_remove = {'#!', 'FIELDS', 'time'}  
    header = [x for x in header if x not in values_to_remove]
    
    df_F = pd.read_csv(filenamesF, sep=r'\s+', skiprows=1, header=None)
    df_F = df_F.drop(df_F.columns[0], axis=1)
    meanF = df_F.mean(axis=0)
    stdvF = df_F.std(axis=0)
    
    df_U = pd.read_csv(filenamesU, sep=r'\s+', skiprows=1, header=None)
    df_U = df_U.drop(df_U.columns[0], axis=1)
    meanU = df_U.mean(axis=0)
    stdvU = df_U.std(axis=0)
    
    # intra- and inter-class covariance
    SwL = stdvF+stdvU
    SwH = 1/((1/stdvF)+(1/stdvU))
    FmU = (meanF-meanU) 
    Sb = (meanF-meanU)*(meanF-meanU)
    lda = Sb/SwL
    hlda = Sb/SwH
    
    df = pd.DataFrame({'labels': header,'meanF': meanF, 'meanU': meanU, 'stdvF': stdvF, 'stdvU': stdvU,
                             'FminusU': FmU, 'Sb': Sb, 'SwH': SwH, 'hlda': hlda, 'SwL': SwL, 'lda': lda})

    return df

# compulsory
parser.add_argument("-s", "--path_plumed", required=True, type=os.path.abspath, help="directory containing the file plumed_final_noduplicate.dat")
parser.add_argument("-f", "--path_folded", required=True, type=os.path.abspath, help="directory containing COLVAR_diff and COLVAR_solvation for the folded trajectory")
parser.add_argument("-u", "--path_unfolded", required=True, type=os.path.abspath, help="directory containing COLVAR_diff and COLVAR_solvation for the unfolded trajectory")
parser.add_argument("-minTs", "--temp_min", required=True, nargs='+', type=float, help="Seven space-separated values for TEMP_MIN, from replica 1 to replica 7. Example --temp_min 319 317 315 310 305 298 290")
parser.add_argument("-maxTs", "--temp_max", required=True, nargs='+', type=float, help="Seven space-separated values for TEMP_MAX, from replica 1 to replica 7. Example --temp_max 325 332 345 360 380 400 420")
parser.add_argument("-p", "--pace", required=True, type=int, help="PACE for the two main CVs (in simulation steps)")
parser.add_argument("-b", "--barrier", required=True, type=float, help="BARRIER for the two main CVs (in kJ/mol)")
parser.add_argument("-rp", "--reference_protein", required=True, type=os.path.abspath, help="reference PDB file for the protein only")
parser.add_argument("-rca", "--reference_CA", required=True, type=os.path.abspath, help="reference PDB file for the CA atoms only")

# optional
parser.add_argument("--outdir", required=False, default="oneopes_files", type=str, help="Name for output directory with OneOPES files. Default: oneopes_files")
parser.add_argument("--pace_minor", required=False, default=100000, type=int, help="PACE for auxiliary CVs (in simulation steps). Default: 100000")
parser.add_argument("--barrier_minor", required=False, default=3, type=float, help="BARRIER for auxiliary CVs (in kJ/mol). Default: 3")
parser.add_argument("--stride", required=False, default=500, type=int, help="STRIDE for COLVAR writing in output. Default: 500")
parser.add_argument("--colvar", required=False, default="COLVAR", type=str, help="Name of the output file containing the collective variables. Default: COLVAR")
parser.add_argument("--opesx_pace", required=False, default=1000, type=int, help="PACE for OPES multithermal bias update. Default: 1000")
parser.add_argument("--opesx_update", required=False, default=0, type=int, help="UPDATE_FROM for starting the OPES multithermal bias application. Default: 0")
parser.add_argument("--opesx_obs", required=False, default=100, type=int, help="OBSERVATION_STEPS for starting the OPES multithermal observation for bias construction, in PACE units. Default: 100")
parser.add_argument("-y", "--yes", required=False, action='store_true', default=False, help='avoid interactivity and run the whole script automatically')

args = parser.parse_args()

args_path_plumed = args.path_plumed
args_path_folded = args.path_folded
args_path_unfolded = args.path_unfolded
args_temp_min = args.temp_min
args_temp_max = args.temp_max
args_pace = args.pace
args_barrier = args.barrier
args_reference_protein = args.reference_protein
args_reference_ca = args.reference_CA
args_outdir = args.outdir
args_pace_minor = args.pace_minor
args_barrier_minor = args.barrier_minor
args_stride = args.stride
args_colvar = args.colvar
args_opesx_pace = args.opesx_pace
args_opesx_update = args.opesx_update
args_opesx_obs = args.opesx_obs
args_yes = args.yes

# Check length of temperatures
if len(args_temp_min) != 7:
    sys.exit(error%('7 values for TEMP_MIN should be given, not ' + str(len(args_temp_min))))
if len(args_temp_max) != 7:
    sys.exit(error%('7 values for TEMP_MAX should be given, not ' + str(len(args_temp_max))))

# minimum sanity check of temperature
for i in range(0, len(args_temp_min)):
    if args_temp_min[i] >= args_temp_max[i]:
        sys.exit(error%(f"The min temperature of replica {i+1} is {args_temp_min[i]}, which is higher than {args_temp_max[i]}, the corresponding max."))

# check that starting plumed file with CVs definition is present
if not os.path.isfile(args_path_plumed + '/plumed_final_noduplicate.dat'):
    sys.exit(error%(f'Missing plumed_final_noduplicate.dat. Please run bioinspired_features.py before running this script or point to the correct directory.'))

# check presence of protein reference files
for file in [args_reference_protein, args_reference_ca]:
    if not os.path.isfile(file):
        sys.exit(error%(f'Could not find {file}.'))

# Check if all COLVAR files needed are present
for file in [[directory,file] for directory in [args_path_folded,args_path_unfolded] for file in ["COLVAR_solvation", "COLVAR_diff"]]:
    if not os.path.isfile(file[0]+"/"+file[1]):
        sys.exit(error%(f'Missing {file[1]} in {file[0]}. Please run bioinspired_features.py before running this script or point to the correct directory.'))

# Check if output directory already exists, if not create it and populate with the things necessary for OneOPES
if os.path.isdir(args_outdir):
    sys.exit(error%(f'Directory {args_outdir} already exists. Delete it or change output directory name with the flag --outdir.'))
os.makedirs(args_outdir)
for i in range(0,8):
    os.makedirs(f'{args_outdir}/rep{i}')
shutil.copyfile(args_reference_protein, f'{args_outdir}/{args_reference_protein.split('/')[-1]}')
shutil.copyfile(args_reference_ca, f'{args_outdir}/{args_reference_ca.split('/')[-1]}')

# some on screen info
print("\n")
print("#########################################################################################")
print("Unbiased folded directory     : " + str(args_path_folded))
print("Unbiased unfolded directory   : " + str(args_path_unfolded))
print("Minimum temperature range     : " + str(args_temp_min))
print("Maximum temperature range     : " + str(args_temp_max))
print("PACE for main CVs             : " + str(args_pace))
print("BARRIER for main CVs          : " + str(args_barrier))
print("Reference file (protein)      : " + str(args_reference_protein))
print("Reference file (protein CA)   : " + str(args_reference_ca))
print("Name of output directory      : " + str(args_outdir))
print("PACE for auxiliary CVs        : " + str(args_pace_minor))
print("BARRIER for auxiliary CVs     : " + str(args_barrier_minor))
print("STRIDE for printing COLVAR    : " + str(args_stride))
print("COLVAR file name for output   : " + str(args_colvar))
print("PACE for OPES multithermal    : " + str(args_opesx_pace))
print("UPDATE_FROM for OPES expanded : " + str(args_opesx_update))
print("OBSERVATION_STEPS for OPES exp: " + str(args_opesx_obs))
print("#########################################################################################")
print("\n")

# check if everything's good with the user
if not args_yes:
    user_input = input("Do you want to continue? (yes/no): ")
    while (user_input.lower() != "yes" and user_input.lower() != "no"):
        print(str(user_input) + " not understood. Please type yes or no.")
        user_input = input("Do you want to continue? (yes/no): ")
    if user_input.lower() == "yes":
        print("Continuing...")
    else:
        print('Exiting ...')
        exit()

# Read COLVAR_diff to extract the main CV
df_diff = read_colvar('COLVAR_diff')
df_diff = df_diff.reset_index()

# Read COLVAR_solvation to extract the auxiliary CV
df_solvation = read_colvar('COLVAR_solvation')
# keeping top 4 carbon atoms
df_solvation_C = df_solvation[df_solvation['labels'].str.startswith('C')]
df_solvation_C = df_solvation_C.sort_values(by=['lda'],ascending=False).head(4).reset_index()
# keeping top 3 carbon atoms
df_solvation_NO = df_solvation[df_solvation['labels'].str.startswith(('N','O'))]
df_solvation_NO = df_solvation_NO.sort_values(by=['lda'],ascending=False).head(3).reset_index()
list_atoms = list(df_solvation_C['labels'])+list(df_solvation_NO['labels'])

for i in range(0,8):
    output = open(f'./{args_outdir}/rep{i}/plumed.dat', 'w')
    for line in open(args_path_plumed + '/plumed_final_noduplicate.dat', "r"):
        if line.lstrip().startswith('PRINT'):
            continue
        if line.lstrip().startswith('MOLINFO'):
            output.write(f"MOLINFO MOLTYPE=protein STRUCTURE=../{args_reference_protein.split('/')[-1]}")
            output.write('\n')
            continue
        if line.lstrip().startswith('rmsd_ca:'):
            output.write(f"rmsd_ca: RMSD TYPE=OPTIMAL REFERENCE=../{args_reference_ca.split('/')[-1]}")
            output.write('\n')
            continue
        output.write(line)
    output.write('\n')
    output.write('# This OneOPES file has been generated by oneopes_files_generation.py on '+ str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'.')
    output.write('\n')
    output.write(f'# Command : python oneopes_files_generation.py -s {args_path_plumed} -f {args_path_folded} -u {args_path_unfolded} -minTs {" ".join(map(str, args_temp_min))} -maxTs {" ".join(map(str, args_temp_max))} -p {args_pace} -b {args_barrier} -rp {args_reference_protein} -rca {args_reference_ca} --outdir {args_outdir} --pace_minor {args_pace_minor} --barrier_minor {args_barrier_minor} --stride {args_stride} --colvar {args_colvar} --opesx_pace {args_opesx_pace} --opesx_update {args_opesx_update} --opesx_obs {args_opesx_obs}')
    output.write('\n')
    output.write('\n')
    output.write('OPES_METAD_EXPLORE ...')
    output.write('\n')
    output.write('\tLABEL=opes')
    output.write('\n')
    output.write('\tARG=' + df_diff.loc[0]['labels'] + ',' + df_diff.loc[1]['labels'])
    output.write('\n')
    output.write('\tSIGMA=' + str(round(df_diff.loc[0]['stdvF'], 2)) + ',' + str(round(df_diff.loc[1]['stdvF'], 2)))
    output.write('\n')
    output.write('\tFILE=Kernels.data')
    output.write('\n')
    output.write('\tSTATE_RFILE=compressed.Kernels')
    output.write('\n')
    output.write('\tSTATE_WFILE=compressed.Kernels')
    output.write('\n')
    output.write('\tPACE=' + str(args_pace))
    output.write('\n')
    output.write('\t#RESTART=YES')
    output.write('\n')
    output.write('\tBARRIER=' + str(args_barrier))
    output.write('\n')
    output.write('... OPES_METAD_EXPLORE')
    if i >= 1:
        output.write('\n')
        output.write('\n')
        output.write('OPES_METAD_EXPLORE ...')
        output.write('\n')
        output.write('\tLABEL=opese1')
        output.write('\n')
        output.write('\tARG=' + str(df_solvation_C.loc[0]['labels'].replace("'", "")))
        output.write('\n')
        output.write('\tSIGMA=' + str(round(df_solvation_C.loc[0]['stdvF'], 2)))
        output.write('\n')
        output.write('\tFILE=Kernelse1.data')
        output.write('\n')
        output.write('\tSTATE_RFILE=compressed.Kernelse1')
        output.write('\n')
        output.write('\tSTATE_WFILE=compressed.Kernelse1')
        output.write('\n')
        output.write('\tPACE=' + str(args_pace_minor))
        output.write('\n')
        output.write('\tBARRIER=' + str(args_barrier_minor))
        output.write('\n')
        output.write('\t#RESTART=YES')
        output.write('\n')
        output.write('... OPES_METAD_EXPLORE')
        output.write('\n')
    if i >= 2:
        output.write('\n')
        output.write('OPES_METAD_EXPLORE ...')
        output.write('\n')
        output.write('\tLABEL=opese2')
        output.write('\n')
        output.write('\tARG=' + str(df_solvation_C.loc[1]['labels'].replace("'", "")))
        output.write('\n')
        output.write('\tSIGMA=' + str(round(df_solvation_C.loc[1]['stdvF'], 2)))
        output.write('\n')
        output.write('\tFILE=Kernelse2.data')
        output.write('\n')
        output.write('\tSTATE_RFILE=compressed.Kernelse2')
        output.write('\n')
        output.write('\tSTATE_WFILE=compressed.Kernelse2')
        output.write('\n')
        output.write('\tPACE=' + str(args_pace_minor))
        output.write('\n')
        output.write('\tBARRIER=' + str(args_barrier_minor))
        output.write('\n')
        output.write('\t#RESTART=YES')
        output.write('\n')
        output.write('... OPES_METAD_EXPLORE')
        output.write('\n')
    if i >= 3:
        output.write('\n')
        output.write('OPES_METAD_EXPLORE ...')
        output.write('\n')
        output.write('\tLABEL=opese3')
        output.write('\n')
        output.write('\tARG=' + str(df_solvation_C.loc[2]['labels'].replace("'", "")))
        output.write('\n')
        output.write('\tSIGMA=' + str(round(df_solvation_C.loc[2]['stdvF'], 2)))
        output.write('\n')
        output.write('\tFILE=Kernelse3.data')
        output.write('\n')
        output.write('\tSTATE_RFILE=compressed.Kernelse3')
        output.write('\n')
        output.write('\tSTATE_WFILE=compressed.Kernelse3')
        output.write('\n')
        output.write('\tPACE=' + str(args_pace_minor))
        output.write('\n')
        output.write('\tBARRIER=' + str(args_barrier_minor))
        output.write('\n')
        output.write('\t#RESTART=YES')
        output.write('\n')
        output.write('... OPES_METAD_EXPLORE')
        output.write('\n')
    if i >= 4:
        output.write('\n')
        output.write('OPES_METAD_EXPLORE ...')
        output.write('\n')
        output.write('\tLABEL=opese4')
        output.write('\n')
        output.write('\tARG=' + str(df_solvation_C.loc[3]['labels'].replace("'", "")))
        output.write('\n')
        output.write('\tSIGMA=' + str(round(df_solvation_C.loc[3]['stdvF'], 2)))
        output.write('\n')
        output.write('\tFILE=Kernelse4.data')
        output.write('\n')
        output.write('\tSTATE_RFILE=compressed.Kernelse4')
        output.write('\n')
        output.write('\tSTATE_WFILE=compressed.Kernelse4')
        output.write('\n')
        output.write('\tPACE=' + str(args_pace_minor))
        output.write('\n')
        output.write('\tBARRIER=' + str(args_barrier_minor))
        output.write('\n')
        output.write('\t#RESTART=YES')
        output.write('\n')
        output.write('... OPES_METAD_EXPLORE')
        output.write('\n')
    if i >= 5:
        output.write('\n')
        output.write('OPES_METAD_EXPLORE ...')
        output.write('\n')
        output.write('\tLABEL=opese5')
        output.write('\n')
        output.write('\tARG=' + str(df_solvation_NO.loc[0]['labels'].replace("'", "")))
        output.write('\n')
        output.write('\tSIGMA=' + str(round(df_solvation_NO.loc[0]['stdvF'], 2)))
        output.write('\n')
        output.write('\tFILE=Kernelse5.data')
        output.write('\n')
        output.write('\tSTATE_RFILE=compressed.Kernelse5')
        output.write('\n')
        output.write('\tSTATE_WFILE=compressed.Kernelse5')
        output.write('\n')
        output.write('\tPACE=' + str(args_pace_minor))
        output.write('\n')
        output.write('\tBARRIER=' + str(args_barrier_minor))
        output.write('\n')
        output.write('\t#RESTART=YES')
        output.write('\n')
        output.write('... OPES_METAD_EXPLORE')
        output.write('\n')
    if i >= 6:
        output.write('\n')
        output.write('OPES_METAD_EXPLORE ...')
        output.write('\n')
        output.write('\tLABEL=opese6')
        output.write('\n')
        output.write('\tARG=' + str(df_solvation_NO.loc[1]['labels'].replace("'", "")))
        output.write('\n')
        output.write('\tSIGMA=' + str(round(df_solvation_NO.loc[1]['stdvF'], 2)))
        output.write('\n')
        output.write('\tFILE=Kernelse6.data')
        output.write('\n')
        output.write('\tSTATE_RFILE=compressed.Kernelse6')
        output.write('\n')
        output.write('\tSTATE_WFILE=compressed.Kernelse6')
        output.write('\n')
        output.write('\tPACE=' + str(args_pace_minor))
        output.write('\n')
        output.write('\tBARRIER=' + str(args_barrier_minor))
        output.write('\n')
        output.write('\t#RESTART=YES')
        output.write('\n')
        output.write('... OPES_METAD_EXPLORE')
        output.write('\n')
    if i >= 7:
        output.write('\n')
        output.write('OPES_METAD_EXPLORE ...')
        output.write('\n')
        output.write('\tLABEL=opese7')
        output.write('\n')
        output.write('\tARG=' + str(df_solvation_NO.loc[2]['labels'].replace("'", "")))
        output.write('\n')
        output.write('\tSIGMA=' + str(round(df_solvation_NO.loc[2]['stdvF'], 2)))
        output.write('\n')
        output.write('\tFILE=Kernelse7.data')
        output.write('\n')
        output.write('\tSTATE_RFILE=compressed.Kernelse7')
        output.write('\n')
        output.write('\tSTATE_WFILE=compressed.Kernelse7')
        output.write('\n')
        output.write('\tPACE=' + str(args_pace_minor))
        output.write('\n')
        output.write('\tBARRIER=' + str(args_barrier_minor))
        output.write('\n')
        output.write('\t#RESTART=YES')
        output.write('\n')
        output.write('... OPES_METAD_EXPLORE')
        output.write('\n')

    ##########################################
    temp_min = args_temp_min
    temp_max = args_temp_max
    ##########################################
    output.write('\n')
    output.write('ene: ENERGY')
    output.write('\n')
    opes_bias = []
    if i >= 1:
        output.write('ecv: ECV_MULTITHERMAL ARG=ene TEMP_MIN=' + str(temp_min[i-1]) + ' TEMP_MAX=' + str(temp_max[i-1]))
        output.write('\n')
        opes_exp_line = f'opesX: OPES_EXPANDED ARG=ecv.* FILE=DeltaFs.data PACE={args_opesx_pace}'
        if args_opesx_update != 0:
            opes_exp_line += f' UPDATE_FROM={args_opesx_update}'
        if args_opesx_obs != 100:
            opes_exp_line += f' OBSERVATION_STEPS={args_opesx_obs}'
        output.write(opes_exp_line)
        output.write('\n')
        output.write('\n')
        output.write('#1 ' + str(temp_min[0]) + '-' + str(temp_max[0]))
        output.write('\n')
        output.write('#2 ' + str(temp_min[1]) + '-' + str(temp_max[1]))
        output.write('\n')
        output.write('#3 ' + str(temp_min[2]) + '-' + str(temp_max[2]))
        output.write('\n')
        output.write('#4 ' + str(temp_min[3]) + '-' + str(temp_max[3]))
        output.write('\n')
        output.write('#5 ' + str(temp_min[4]) + '-' + str(temp_max[4]))
        output.write('\n')
        output.write('#6 ' + str(temp_min[5]) + '-' + str(temp_max[5]))
        output.write('\n')
        output.write('#7 ' + str(temp_min[6]) + '-' + str(temp_max[6]))

        opes_bias = ['opesX.bias']
    output.write('\n')
    output.write('\n')
    for j in range(1,i+1):
        temp = 'opese'+ str(j) + '.bias'
        opes_bias.append(temp)
    print_statement = f'PRINT STRIDE={args_stride} ARG=opes.bias,ene,rmsd_ca,' + df_diff.loc[0]['labels'] + ',' + df_diff.loc[1]['labels']
    if i > 0:
        print_statement = print_statement + ',' + str(list(list_atoms)) + ',' + str(list(opes_bias)) 
    print_statement += f' FILE={args_colvar}'
    print_statement = print_statement.replace("]"," ").replace("["," ").replace(", ",",").replace(" ,",",").replace("= ","=").replace("'","")
    output.write(print_statement)
    output.write('\n')
    output.write('\n') 
    output.close()
print("Done !")
