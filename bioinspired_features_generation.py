import argparse
from argparse import ArgumentParser
from pathlib import Path
import mdtraj
import numpy as np
import pandas as pd
import re
import os
import sys
import subprocess
import logging
import datetime

version = 1.0
error = '--- ERROR: %s \n'
logging.basicConfig(filename='bioinspired_features.log', filemode='a', format='%(levelname)s:%(message)s', level=logging.INFO)

# nomenclatures for hbond donors and acceptors
# tested agains desamber, amber99SB-ildn, charmm27, GROMOS, OPLS
# definition of acceptors wrt mdtraj nomenclature
atom_hacceptor_dict= {
    "O": "C",
    "OD1": "CG",
    "OD2": "CG",
    "OE1": "CD",
    "OE2": "CD",
    "OG": "CB",
    "OG1": "CB",
    "SD": ["CG","CE"],
    "SG": "CB",
    "OH": "CZ",
    "OXT": "C",
    "ND1": ["CG", "CE1"],
    "NE2": ["CD2", "CE1"]
    }
hbond_acceptor = list(atom_hacceptor_dict)

# definition of donors wrt mdtraj nomenclature
hatom_hdonor_dict = {
    "H": "N",
    "H2": "N",
    "H3": "N",
    "HD1": "ND1",
    "HE1": "NE1",
    "HE2": "NE2",
    "HZ1": "NZ",
    "HZ2": "NZ",
    "HZ3": "NZ",
    "HE21": "NE2",
    "HE22": "NE2",
    "HE": "NE",
    "HG": "OG",
    "HG": "SG",
    "HH11": "NH1",
    "HH12": "NH1",
    "HH21": "NH2",
    "HH22": "NH2",
    "HD21": "ND2",
    "HD22": "ND2",
    "HG1": "OG1",
    "HH": "OH"
    }
hbond_donor = list(hatom_hdonor_dict)

# Function definitions
def read_colvar(name: str, cv_prefix: str):
    """
    Read the colvar files of the folded and unfolded directories
    Returns a dataframe with some statistics on loaded CVs (except time)
    """
    # Reads COLVAR file
    filenamesF = args_folded_dir + "/" + name
    filenamesU = args_unfolded_dir + "/" + name

    header = pd.read_csv(filenamesF, sep=r'\s+', nrows=0).columns.tolist()
    values_to_remove = {'#!', 'FIELDS', 'time'}
    header = [x for x in header if x not in values_to_remove]

    df_F = pd.read_csv(filenamesF, sep=r'\s+',skiprows=1, header=None)
    df_F = df_F.drop(df_F.columns[0], axis=1)
    meanF = df_F.mean(axis=0)
    stdvF = df_F.std(axis=0)

    df_U = pd.read_csv(filenamesU, sep=r'\s+',skiprows=1, header=None)
    df_U = df_U.drop(df_U.columns[0], axis=1)
    meanU = df_U.mean(axis=0)
    stdvU = df_U.std(axis=0)

    # intra- and inter-class covariance
    SwL = stdvF + stdvU
    SwH = 1/((1/stdvF)+(1/stdvU))
    FmU = (meanF-meanU)
    Sb = (meanF-meanU)*(meanF-meanU)
    lda = Sb/SwL
    hlda = Sb/SwH

    df = pd.DataFrame({'labels': header,'meanF': meanF, 'meanU': meanU, 'stdvF': stdvF, 'stdvU': stdvU,
                             'FminusU': FmU, 'Sb': Sb, 'SwH': SwH, 'hlda': hlda, 'SwL': SwL, 'lda': lda})

    df = df[df['labels'].str.startswith(cv_prefix)]
    return df


# Parse user inputs
parser = argparse.ArgumentParser(
    prog = "bioinspired_features",
    description = "Designs and filters bioinspired features for use in enhanced sampling simulations of peptide folding.",
    epilog = "Thanks for using %(prog)s! ",
)

# compulsory
parser.add_argument("-F", "--folded", required=True, type=os.path.abspath, help="unbiased folded trajectory (xtc file)")
parser.add_argument("-U", "--unfolded", required=True, type=os.path.abspath, help="unbiased unfolded trajectory (xtc file)")
parser.add_argument("-r", "--reference_pdb_solvated", required=True, type=os.path.abspath, help="reference PDB file of the solvated system")
parser.add_argument("-rp", "--reference_protein", required=True, type=os.path.abspath, help="reference PDB file for the protein only")
parser.add_argument("-rca", "--reference_CA", required=True, type=os.path.abspath, help="reference PDB file for the CA atoms only")
parser.add_argument("-mc", "--mcfile", required=True, type=os.path.abspath, help="PLUMED mcfile")

# other options
parser.add_argument("-l", "--lda", required=False, default=0.3, type=float, help="lda value to use for the filtering. Default: 0.3")
parser.add_argument("-c", "--cutoff", required=False, default=0.6, type=float, help="cutoff value to use for the filtering. Default: 0.6nm")
parser.add_argument("-s", "--stride", required=False, default=10, type=int, help="STRIDE for the PLUMED file in the filtering steps. Default: 10")
parser.add_argument("-e", "--explicit", required=False, action='store_true', default=False, help='writes features in the explicit fashion')
parser.add_argument("-py", "--pymol", required=False, action='store_true', default=False, help='saves a session of pymol with the main hydrogen bonds highlighted; requires the pymol module')
parser.add_argument("-y", "--yes", required=False, action='store_true', default=False, help='avoid interactivity and run the whole script automatically')

args = parser.parse_args()

inp_dir = os.getcwd()
# folded
args_path_folded = args.folded
args_folded_dir = os.path.dirname(args_path_folded)
args_folded_trajectory = os.path.basename(args_path_folded)
# unfolded
args_path_unfolded = args.unfolded
args_unfolded_dir = os.path.dirname(args_path_unfolded)
args_unfolded_trajectory = os.path.basename(args_path_unfolded)
# references
args_reference = args.reference_pdb_solvated
args_reference_protein = args.reference_protein
args_reference_ca = args.reference_CA
args_mcfile = args.mcfile
# additional parameters
args_lda = args.lda
args_cutoff = args.cutoff
args_stride = args.stride
args_explicit = args.explicit
args_pymol = args.pymol
args_yes = args.yes

# Check format of input
# Should we check also the mcfile?
if (args_folded_trajectory[-3:] != "xtc" or args_unfolded_trajectory[-3:] != "xtc"):
    sys.exit(error%('Trajectories must be an xtc file.'))
if (args_reference[-3:] != "pdb" or args_reference_protein[-3:] != "pdb" or args_reference_ca[-3:] != "pdb"):
    sys.exit(error%('Reference structures must be a pdb file.'))
# Check files availability
missing_files = []
for infile in [args_reference_protein, args_reference, args_reference_ca, args_mcfile, args_path_folded, args_path_unfolded]:
    if not os.path.isfile(infile):
        missing_files.append(infile)
if len(missing_files) > 0:
    sys.exit(error%('The following input files have not been found:\n\n  > ' + '\n  > '.join(missing_files)))
# check plumed sourcing
try:
    subprocess.run(["plumed", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
except:
    sys.exit(error%('PLUMED is not installed or has not been sourced properly.'))

# load reference PDB file
topology_solv = mdtraj.load(args_reference).topology
ref_solv, bonds_solv = topology_solv.to_dataframe()

# checks the type of water (# atoms in the model)
# safer to check the number of atoms rather than the names as we pick up also broken pdbs
n_atoms_water = len(topology_solv.select("water"))
n_res_waters = 0
for residue in topology_solv.residues:
    if residue.is_water: n_res_waters += 1
water_model = n_atoms_water/n_res_waters

# checks that the number is divisible by three or four, not perfect but relatively robust
if (water_model == 3.0 or water_model == 4.0): # 3 and 4 points supported
    water_model = int(water_model)
elif (water_model % 1 == 0):
    sys.exit(error%('The water model detected has ' + str(water_model) + ' atoms per residue, which is not supported'))
else:
    sys.exit(error%('Non-integer number of water atoms detected (' + str(water_model) + ' atoms per residue). Check the reference file: ' + args_reference))

# check if pymol works if the saving of a pymol session is requested
if args_pymol:
    try:
        import pymol
    except ImportError:
        print("WARNING: Module pymol not found. Is it installed?")
        print("WARNING: Will not print the pymol session...")
        pymol_viz = False
    except:
        print("WARNING: Cannot use pymol module, but it was found. Something else went wrong!")
        print("WARNING: Will not print the pymol session...")
        pymol_viz = False
    else:
        from pymol import cmd

# store first and last water oxygen atoms
first_oxygen = list(ref_solv[(ref_solv['resName'] == 'HOH') & (ref_solv['name'] == 'O')].head(1)['serial'])[0]
last_oxygen = list(ref_solv[(ref_solv['resName'] == 'HOH') & (ref_solv['name'] == 'O')].tail(1)['serial'])[0]

# some on screen info
print("\n")
print("################################################################################################################")
print("Working directory           : " + str(inp_dir))
print("Reference file (solvated)   : " + str(args_reference))
print("Reference file (protein)    : " + str(args_reference_protein))
print("Reference file (protein CA) : " + str(args_reference_ca))
print("Folded trajectory file      : " + str(args_path_folded))
print("Unfolded trajectory file    : " + str(args_path_unfolded))
print("Reference mcfile            : " + str(args_mcfile))
print("Water model detected        : " + str(water_model) + "-points water model.")
print("Stride                      : " + str(int(args_stride)))
print("LDA cutoff value            : " + str(args_lda))
print("Cutoff value for the H-bonds: " + str(args_cutoff))
print("Writing explicit features   : " + str(args_explicit))
print("Pymol session generation    : " + str(args_pymol))
if args_pymol: print(f"Pymol session saved in file : {inp_dir}/symmary_pymol_session.pse")
print("################################################################################################################")
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

logging.info(' BIOINSPIRED_FEATURES_GENERATION Python script version ' + str(version))
logging.info(' Please read and cite the following reference:')
logging.info(' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx')
logging.info(' Main contributors: Margaux Héritier, Valerio Rizzi, Nicola Piasentin, Simone Aureli, Francesco Luigi Gervasio.')
logging.info(' Command Line:')
logging.info(' python bioinspired_features_generation_v.1.py -F '+ str(args_path_folded) + ' -U ' + str(args_path_unfolded) + ' -r ' + str(args_reference) + ' -rp ' + str(args_reference_protein) + ' -rca ' + str(args_reference_ca) + ' -mc ' + str(args_mcfile) + ' -l ' + str(args_lda) + ' -c ' + str(args_cutoff) + ' -s ' + str(args_stride) + (' -e ' if args_explicit else '') + (' -py ' if args_pymol else '')  + (' -y ' if args_yes else ''))
logging.info(' Started on '+ str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'.')
logging.info('\n')

# load reference PDB file (only protein)
topology = mdtraj.load(args_reference_protein).topology
table, bonds = topology.to_dataframe()

# generate a list of dataframes of all the hydrogens that could participate in H-bonds in args_reference_protein
list_hbond_donors = []
for atom in hbond_donor:
    temp = table[table['name'] == atom].reset_index()
    if len(temp) != 0: list_hbond_donors.append(temp)

# generate a list of dataframes of all the H-bond acceptors that could participate in H-bonds in args_reference_protein
list_hbond_acceptors = []
for atom in hbond_acceptor:
    temp = table[table['name'] == atom].reset_index()
    if len(temp) != 0: list_hbond_acceptors.append(temp)

###################################################################################
# First filter to keep only the relevant contacts

# write a PLUMED script to get the distances of all possible H-bonds in the structure
print("Writing plumed_1.dat file ...")
output= open('plumed_1.dat', 'w')
output.write(f"MOLINFO MOLTYPE=protein STRUCTURE={args_reference_protein}")
output.write('\n')
output.write('# all distances between possible H-bonds')
output.write('\n')

# print out all the combination of possible intra-protein hbonds
for i in range(0, len(list_hbond_donors)): # to go through the H-bond donors list
    for j in range(0, list_hbond_donors[i].shape[0]): 
        for z in range(0, len(list_hbond_acceptors)):
            for y in range(0, list_hbond_acceptors[z].shape[0]): # to go through the H-bond acceptors list
                label = str(list_hbond_donors[i]['name'][j]) + "_" + str(list_hbond_donors[i]['serial'][j]) + "-" + str(list_hbond_acceptors[z]['name'][y]) + "_" + str(list_hbond_acceptors[z]['serial'][y])
                s = label + ": DISTANCE ATOMS=" + str(list_hbond_donors[i]['serial'][j]) + "," + str(list_hbond_acceptors[z]['serial'][y])
                output.write(s)
                output.write('\n')
            
output.write('\n')
# compute the distance each args_stride frames to save computational time since there will be a high number of distances to compute
output.write('PRINT STRIDE=' + str(args_stride) + ' ARG=* FILE=COLVAR_first_filter') 
output.write('\n')
output.close()

print("Running plumed_1.dat on the folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_1.dat --mf_xtc {args_folded_trajectory} --pdb {args_reference} 1> plumed_1_folded.out")
print("Running plumed_1.dat on the unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_1.dat --mf_xtc {args_unfolded_trajectory} --pdb {args_reference} 1> plumed_1_unfolded.out")
os.chdir(inp_dir)

print("Reading COLVAR file for the first filter ...")
df = read_colvar('COLVAR_first_filter', '')
# drop H bonds with average length < 0.4nm
df = df.drop(df[(df['meanF'] > 0.4) & (df['meanU'] > 0.4) ].index).dropna()
# to remove the potential contacts within the same residues that don't display any difference between folded and unfolded
# drop H bonds with average difference between folded and unfolded < 0.01 nm
df = df.drop(df[(df['FminusU'].abs() < 0.01)].index).dropna()
df = df.sort_values(by=['labels'], ascending=True).reset_index()

###############################################################################################3
# Second filter to keep only the relevant contacts

# writes PLUMED script to obtain the angle information from the relevant H-bonds in df
print("Writing plumed_2.dat file to obtain angle information ...")
output = open('plumed_2.dat', 'w')
output.write(f"MOLINFO MOLTYPE=protein STRUCTURE={args_reference_protein}")
output.write('\n')
output.write('\n')

for i in range(0,len(df)):
    labels = df['labels'][i]
    match_atomtype = re.findall(r'([A-Za-z0-9]+)_\d+', labels)
    match_number = re.findall(r'_(\d+)', labels)

    A = match_atomtype[1] # the atom type of A
    A_serial = match_number[1] # the atom number of A
    
    H = match_atomtype[0] # the atom type of H
    H_serial = match_number[0] # the atom number of H

    D = hatom_hdonor_dict[match_atomtype[0]] # the atom type of D, corresponding to the match in the dictionnary
    D_resSeq_temp = table[table['serial'] == int(match_number[0])]['resSeq'].reset_index() # temporary variable, dataframe of the atom number of H
    D_resSeq = D_resSeq_temp['resSeq'].loc[0] # residue number of D, which is the same one as H
    D_serial_temp = table[(table['resSeq'] == int(D_resSeq)) & (table['name'] == D)]['serial'].reset_index() # temporary variable, dataframe of the residue corresponding to the residue number of D and its atom type
    if D_serial_temp.empty: # if temp is empty, it means that the hydrogen atom cannot form an H-bond because it is not covalently linked to an atomtype in the dictionnary
        continue
    else:  
        D_serial = D_serial_temp['serial'].loc[0] # atom number of D
    
    B = atom_hacceptor_dict[match_atomtype[1]]
    
    if isinstance(B,str): # meaning it is a string so B is only one value
        B_resSeq_temp = table[table['serial'] == int(match_number[1])]['resSeq'].reset_index()
        B_resSeq = B_resSeq_temp['resSeq'].loc[0]
        B_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B)]['serial'].reset_index()
        B_serial = B_serial_temp['serial'].loc[0]

        # creating the labels that will be used in PLUMED
        label_cont = "cont_"+str(df['labels'][i])
        label_HD = "H-D_"+str(H)+"_"+str(H_serial)+"-"+str(D)+"_"+str(D_serial)
        label_HA = "H-A_"+str(H)+"_"+str(H_serial)+"-"+str(A)+"_"+str(A_serial)
        label_HB = "H-B_"+str(H)+"_"+str(H_serial)+"-"+str(B)+"_"+str(B_serial)
        label_DA = "D-A_"+str(D)+"_"+str(D_serial)+"-"+str(A)+"_"+str(A_serial)
        label_AB = "A-B_"+str(A)+"_"+str(A_serial)+"-"+str(B)+"_"+str(B_serial)
        label_ang_DHA = "ang_DHA_"+str(labels)
        label_ang_BAH = "ang_BAH_"+str(labels)
        label_hbond = "hbond_"+str(labels)

        # writing the different ingredients that make the H-bond
        s_cont=label_cont+": COORDINATION GROUPA="+str(match_number[0])+" GROUPB="+str(match_number[1])+" SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=8}"
        sHD=label_HD+": DISTANCE ATOMS="+str(H_serial)+","+str(D_serial)
        sHA=label_HA+": DISTANCE ATOMS="+str(H_serial)+","+str(A_serial)
        sHB=label_HB+": DISTANCE ATOMS="+str(H_serial)+","+str(B_serial)
        sDA=label_DA+": DISTANCE ATOMS="+str(D_serial)+","+str(A_serial)
        sAB=label_AB+": DISTANCE ATOMS="+str(A_serial)+","+str(B_serial)
        sang_DHA = label_ang_DHA+": CUSTOM ARG="+label_DA+","+label_HD+","+label_HA+" FUNC=x/(y+z) PERIODIC=NO"
        sang_BAH = label_ang_BAH+": CUSTOM ARG="+label_HB+","+label_HA+","+label_AB+" FUNC=x/(y+z) PERIODIC=NO"
        shbond=label_hbond+": CUSTOM ARG="+label_cont+","+label_ang_DHA+","+label_ang_BAH+" FUNC=x*y*z PERIODIC=NO"
        
        output.write("# Distance between pairs of atoms of residues "+str(D_resSeq)+" and "+str(B_resSeq))
        output.write('\n')
        output.write(s_cont)
        output.write('\n')
        output.write(sHD)
        output.write('\n')
        output.write(sHA)
        output.write('\n')
        output.write(sHB)
        output.write('\n')
        output.write(sDA)
        output.write('\n')
        output.write(sAB)
        output.write('\n')
        output.write(sang_DHA)
        output.write('\n')
        output.write(sang_BAH)
        output.write('\n')
        output.write(shbond)
        output.write('\n')
        output.write('\n')

    else: # meaning it is ND1 from histidine or SD from methionine that are bound to two different B => COM between two B
        B_resSeq_temp = table[table['serial'] == int(match_number[1])]['resSeq'].reset_index()
        B_resSeq = B_resSeq_temp['resSeq'].loc[0]

        B_1 = B[0]
        B_2 = B[1]
        
        B_1_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B_1)]['serial'].reset_index()
        B_1_serial = B_1_serial_temp['serial'].loc[0]

        B_2_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B_2)]['serial'].reset_index()
        B_2_serial = B_2_serial_temp['serial'].loc[0]

        # creating the labels that will be used in PLUMED
        label_cont = "cont_"+str(df['labels'][i])
        label_com = "com_"+str(B_1)+"_"+str(B_1_serial)+"-"+str(B_2)+"_"+str(B_2_serial)
        label_HD = "H-D_"+str(H)+"_"+str(H_serial)+"-"+str(D)+"_"+str(D_serial)
        label_HA = "H-A_"+str(H)+"_"+str(H_serial)+"-"+str(A)+"_"+str(A_serial)
        label_HB = "H-B_"+str(H)+"_"+str(H_serial)+"-"+str(B_1)+"_"+str(B_1_serial) # for the sake of the nedt steps I kept only B_1 for the name
        label_DA = "D-A_"+str(D)+"_"+str(D_serial)+"-"+str(A)+"_"+str(A_serial)
        label_AB = "A-B_"+str(A)+"_"+str(A_serial)+"-"+str(B_1)+"_"+str(B_1_serial)
        label_ang_DHA = "ang_DHA_"+str(labels)
        label_ang_BAH = "ang_BAH_"+str(labels)
        label_hbond = "hbond_"+str(labels)

        # writing the different ingredients that make the H-bond
        s_cont=label_cont+": COORDINATION GROUPA="+str(match_number[0])+" GROUPB="+str(match_number[1])+" SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=8}"
        com = label_com+": COM ATOMS="+str(B_1_serial)+","+str(B_2_serial)
        sHD=label_HD+": DISTANCE ATOMS="+str(H_serial)+","+str(D_serial)
        sHA=label_HA+": DISTANCE ATOMS="+str(H_serial)+","+str(A_serial)
        sHB=label_HB+": DISTANCE ATOMS="+str(H_serial)+","+str(label_com)
        sDA=label_DA+": DISTANCE ATOMS="+str(D_serial)+","+str(A_serial)
        sAB=label_AB+": DISTANCE ATOMS="+str(A_serial)+","+str(label_com)
        sang_DHA = label_ang_DHA+": CUSTOM ARG="+label_DA+","+label_HD+","+label_HA+" FUNC=x/(y+z) PERIODIC=NO"
        sang_BAH = label_ang_BAH+": CUSTOM ARG="+label_HB+","+label_HA+","+label_AB+" FUNC=x/(y+z) PERIODIC=NO"
        shbond=label_hbond+": CUSTOM ARG="+label_cont+","+label_ang_DHA+","+label_ang_BAH+" FUNC=x*y*z PERIODIC=NO"
        
        output.write("# Distance between pairs of atoms of residues "+str(D_resSeq)+" and "+str(B_resSeq))
        output.write('\n')
        output.write(s_cont)
        output.write('\n')
        output.write(com)
        output.write('\n')
        output.write(sHD)
        output.write('\n')
        output.write(sHA)
        output.write('\n')
        output.write(sHB)
        output.write('\n')
        output.write(sDA)
        output.write('\n')
        output.write(sAB)
        output.write('\n')
        output.write(sang_DHA)
        output.write('\n')
        output.write(sang_BAH)
        output.write('\n')
        output.write(shbond)
        output.write('\n')
        output.write('\n')
    
        output.write('\n')

output.write('\n')
output.write('\n')
output.write('PRINT STRIDE='+str(args_stride)+ ' ARG=* FILE=COLVAR_second_filter')  
output.close()

# if the same hydrogen participates in different H-bonds, then the definition to its distance H-D will be printed multiple times. Therefore we need to eliminate the duplicates.
lines_seen = set() # holds lines already seen
outfile = open('plumed_2_noduplicate.dat', "w")
for line in open('plumed_2.dat', "r"):
    if line not in lines_seen: # not a duplicate
        outfile.write(line)
        lines_seen.add(line)
outfile.close()

print("Running plumed_2.dat file on folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_2_noduplicate.dat --mf_xtc {args_folded_trajectory} --pdb {args_reference} 1> plumed_2_folded.out")
print("Running plumed_2.dat file on unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_2_noduplicate.dat --mf_xtc {args_unfolded_trajectory} --pdb {args_reference} 1> plumed_2_unfolded.out")
os.chdir(inp_dir)

print("Reading COLVAR files ...")
df2 = read_colvar('COLVAR_second_filter','hbond_')
# df2.sort_values(by=['meanF'],ascending=False)
# extract H-bonds that are significant in folded state
df2['Significant folded'] = (df2['FminusU'] > 0) & (df2['lda'] > args_lda)
# edtract H-bonds that are significant in unfolded state
df2['Significant unfolded'] = (df2['FminusU'] < 0) & (df2['lda'] > args_lda)
# remove contacts that are not significant in neither 
df2 = df2.drop(df2[(df2['Significant folded'] == False) & (df2['Significant unfolded'] == False) ].index).dropna()
# define the hard H-bonds
df2['Hard'] = ((df2['Significant folded'] == True) & (df2['meanF'] >= args_cutoff)) | ((df2['Significant unfolded'] == True) & (df2['meanU'] >= args_cutoff))
df_hard_sigF_sigU = df2[df2['Hard']== True]
df_hard_sigF_sigU = df_hard_sigF_sigU.sort_values(by=['lda'],ascending=False).reset_index(drop=True)
# create dataframes for the hard H-bonds depending on their significance
df_hard_sigF = df_hard_sigF_sigU[df_hard_sigF_sigU['Significant folded']== True]
df_hard_sigF = df_hard_sigF.sort_values(by=['lda'],ascending=False).reset_index(drop=True)
df_hard_sigU = df_hard_sigF_sigU[df_hard_sigF_sigU['Significant unfolded']== True]
df_hard_sigU = df_hard_sigU.sort_values(by=['lda'],ascending=False).reset_index(drop=True)

# to obtain the soft H-bonds
df3 = read_colvar('COLVAR_second_filter','cont_')
# remove the "hard" H-bonds from this new dataframe
labels_to_drop = []

for i in range(len(df_hard_sigF_sigU)):
    labels = df_hard_sigF_sigU['labels'][i]
    label_clean = labels.replace('hbond_', '')
    temp = 'cont_' + label_clean
    labels_to_drop.append(temp)

df_soft = df3[~df3['labels'].isin(labels_to_drop)].copy()
df_soft.loc[:, 'labels'] = df_soft['labels'].apply(lambda d: d.replace('cont_', 'hbond_'))

# create dataframes for the soft H-bonds depending on their significance
df_soft.loc[:, 'Significant folded'] = (df_soft['FminusU'] > 0) & (df_soft['meanF'] > args_cutoff) & (df_soft['lda'] > args_lda)
df_soft.loc[:, 'Significant unfolded'] = (df_soft['FminusU'] < 0) & (df_soft['meanU'] > args_cutoff) & (df_soft['lda'] > args_lda)
df_soft_sigF = df_soft[df_soft['Significant folded']== True]
df_soft_sigF = df_soft_sigF.sort_values(by=['lda'],ascending=False).reset_index(drop=True)
df_soft_sigU = df_soft[df_soft['Significant unfolded']== True]
df_soft_sigU = df_soft_sigU.sort_values(by=['lda'],ascending=False).reset_index(drop=True)
df_soft_sigF_sigU = pd.concat([df_soft_sigF,df_soft_sigU]).reset_index(drop=True)

# concatenating a dataframe with all the hard and soft H-bonds
df_hard_soft = pd.concat([df_hard_sigF_sigU,df_soft_sigF_sigU]).reset_index(drop=True)
# fill in the 'Hard' label with false for soft hbonds
df_hard_soft['Hard'] = df_hard_soft['Hard'].notna()

#########################################################################################3
# Generating the exclusion list for the H-bonds

# script to identify which hydrogen and H-bond acceptor should be excluded.
# This script will compute the distance between a virtual atom located 2.5 A from either H or A .
# We're not going to use the virtual atoms directly but this is a necessary step to come up with the exclusion list

# generates dataframe of all hydrogens and acceptors from folded.pdb
all_H = table[table['element'] == "H"].reset_index()
all_acceptors = table[(table['element'] == "O") | (table['element'] == "N")].reset_index()

print("Writing plumed_3.dat ...")
output = open('plumed_3.dat', 'w')
output.write(f"MOLINFO MOLTYPE=protein STRUCTURE={args_reference_protein}")
output.write('\n')

for i in range(0,len(df_hard_soft)):
    labels = df_hard_soft['labels'][i]
    match_atomtype = re.findall(r'([A-Za-z0-9]+)_\d+', labels)
    match_number = re.findall(r'_(\d+)', labels)
    original_label = labels.replace("hbond_", "")
    original_label = original_label.replace("cont_", "")

    A = match_atomtype[1] # the atom type of A
    A_serial = match_number[1] # the atom number of A

    H = match_atomtype[0] # the atom type of H
    H_serial = match_number[0] # the atom number of H
    
    D = hatom_hdonor_dict[match_atomtype[0]] # the atom type of D, corresponding to the match in the dictionnary
    D_resSeq_temp = table[table['serial'] == int(match_number[0])]['resSeq'].reset_index() # temporary variable, dataframe of the atom number of H
    D_resSeq = D_resSeq_temp['resSeq'].loc[0] # residue number of D, which is the same one as H
    D_serial_temp = table[(table['resSeq'] == int(D_resSeq)) & (table['name'] == D)]['serial'].reset_index() # temporary variable, dataframe of the residue corresponding to the residue number of D and its atom type
    D_serial = D_serial_temp['serial'].loc[0] # atom number of D
    
    B = atom_hacceptor_dict[match_atomtype[1]]

    if isinstance(B,str): # meaning it is a string so B is only one value
        B_resSeq_temp = table[table['serial'] == int(match_number[1])]['resSeq'].reset_index()
        B_resSeq = B_resSeq_temp['resSeq'].loc[0]
        B_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B)]['serial'].reset_index()
        B_serial = B_serial_temp['serial'].loc[0]

        label_VA = "VA_"+str(original_label)
        label_VH = "VH_"+str(original_label)
    
        s_vA = label_VA +": GHOST ATOMS="+str(B_serial)+","+str(A_serial)+","+str(int(A_serial)+1)+" COORDINATES=0.25,0.0,0.0"  
        s_VH = label_VH +": GHOST ATOMS="+str(D_serial)+","+str(H_serial)+","+str(int(H_serial)+1)+" COORDINATES=0.25,0.0,0.0"
    
        output.write('\n')
        output.write(s_vA)
        output.write('\n')
        output.write(s_VH)
        output.write('\n')

    
        for i in range(0,len(all_acceptors)): # for hbond_donor list   ##DEBUG, switching all_H with all_acceptors
                label_ghost_donor = "virtual_HD_"+str(original_label)+"-"+str(all_acceptors['name'][i])+"_"+str(all_acceptors['serial'][i])
                s_ghost_donor=label_ghost_donor+": DISTANCE ATOMS="+label_VH+","+str(all_acceptors['serial'][i])
                output.write(s_ghost_donor)
                output.write('\n')
    
        output.write('\n')
    
        for d in range(0,len(all_H)):          ##DEBUG 
                label_ghost_acceptor = "virtual_AB_"+str(original_label)+"-"+str(all_H['name'][d])+"_"+str(all_H['serial'][d])
                s_ghost_acceptor=label_ghost_acceptor+": DISTANCE ATOMS="+label_VA+","+str(all_H['serial'][d])
                output.write(s_ghost_acceptor)
                output.write('\n')
        
        output.write('\n')

    else: # meaning it is ND1 from histidine or SD from methionine that are bound to two different B => COM between two B
        B_resSeq_temp = table[table['serial'] == int(match_number[1])]['resSeq'].reset_index()
        B_resSeq = B_resSeq_temp['resSeq'].loc[0]
    
        B_1 = B[0]
        B_2 = B[1]
            
        B_1_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B_1)]['serial'].reset_index()
        B_1_serial = B_1_serial_temp['serial'].loc[0]
    
        B_2_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B_2)]['serial'].reset_index()
        B_2_serial = B_2_serial_temp['serial'].loc[0]

        label_com = "com_"+str(B_1)+"_"+str(B_1_serial)+"-"+str(B_2)+"_"+str(B_2_serial)
        com = label_com+": COM ATOMS="+str(B_1_serial)+","+str(B_2_serial)

        label_VA = "VA_"+str(original_label)
        label_VH = "VH_"+str(original_label)
        
        s_vA = label_VA +": GHOST ATOMS="+str(com)+","+str(A_serial)+","+str(int(A_serial)+1)+" COORDINATES=0.25,0.0,0.0"  
        s_VH = label_VH +": GHOST ATOMS="+str(D_serial)+","+str(H_serial)+","+str(int(H_serial)+1)+" COORDINATES=0.25,0.0,0.0"
        
        output.write('\n')
        output.write(s_vA)
        output.write('\n')
        output.write(s_VH)
        output.write('\n')
    
        for i in range(0,len(all_acceptors)): # for hbond_donor list  ##DEBUG
            label_ghost_donor = "virtual_HD_"+str(original_label)+"-"+str(all_acceptors['name'][i])+"_"+str(all_acceptors['serial'][i])
            s_ghost_donor=label_ghost_donor+": DISTANCE ATOMS="+label_VH+","+str(all_acceptors['serial'][i])
            output.write(s_ghost_donor)
            output.write('\n')
        
        output.write('\n')
        
        for d in range(0,len(all_H)):   ##DEBUG
            label_ghost_acceptor = "virtual_AB_"+str(original_label)+"-"+str(all_H['name'][d])+"_"+str(all_H['serial'][d])
            s_ghost_acceptor=label_ghost_acceptor+": DISTANCE ATOMS="+label_VA+","+str(all_H['serial'][d])
            output.write(s_ghost_acceptor)
            output.write('\n')
            
        output.write('\n')
        output.write('\n')
    
output.write('\n')
output.write('\n')
output.write('PRINT STRIDE='+str(args_stride)+ ' ARG=* FILE=COLVAR_third_filter')  
output.write('\n')
output.close()

print("Running plumed_3.dat on folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_3.dat --mf_xtc {args_folded_trajectory}  --pdb {args_reference} 1> plumed_3_folded.out")
print("Running plumed_3.dat on unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_3.dat --mf_xtc {args_unfolded_trajectory}  --pdb {args_reference} 1> plumed_3_unfolded.out")
os.chdir(inp_dir)

# read the COLVAR files for the virtual contacts
df_virtual = read_colvar('COLVAR_third_filter','virtual_')
df_virtual = df_virtual.sort_values(by=['meanF'],ascending=True)
# here we are keeping only the contacts that are more than 4 A away from the H-bond in the folded state
df_virtual = df_virtual[df_virtual['meanF'] > 0.4].dropna() 
# df_virtual.sort_values(by=['meanF'],ascending=True)

#############################################################################
# Adding the side chains

# writes a PLUMED file for the side-chain contacts
print("Writing plumed_SC.dat for the side-chain contacts ...")
output = open('plumed_SC.dat', 'w')
output.write(f"MOLINFO MOLTYPE=protein STRUCTURE={args_reference_protein}")
output.write('\n')
output.write('\n')

# to save the amount of residues forming the protein
target = mdtraj.load(args_reference_protein)
numres=target.topology.select("name CA").shape

# to calculate the geometric center of the side chains (SC)
resID_nogly=[]
for i in range(0,numres[0]):
    p=target.topology.select("resid "+str(i)+" and not backbone and not (type H)")+1
    p = p.tolist()
    if len(p) !=0:
        s="SC"+str(i+1)+": CENTER ATOMS="+str(p)+"MASS"
        output.write(s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","="))
        output.write('\n')
    else:
        p=target.topology.select("resid "+str(i)+" and backbone and name CA")+1
        s="SC"+str(i+1)+": CENTER ATOMS="+str(p)+"MASS"
        output.write(s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","="))
        output.write('\n')
output.write('\n')

# CV of the distances between the SCs (adjacent SC are not taken into account)
# the list resID_nogly is needed to take into account the removal of GLYs
for i in range(1,numres[0]+1):
    for j in range(i+2,numres[0]+1):
        label = "contside"+str(i)+"-"+str(j)+": "
        output.write(str(label)+"COORDINATION GROUPA=SC"+str(i)+" GROUPB=SC"+str(j)+"  SWITCH={RATIONAL D_0=0.0 R_0=0.80 NN=4 MM=8}")
        output.write('\n')
output.write('\n')

# calculate the distance between the two center of masse of residues pairs
for i in range(1,numres[0]+1):
    for j in range(i+2,numres[0]+1):
        label = "SC"+str(i)+"-"+str(j)+": "
        output.write(str(label)+"CENTER ATOMS=SC"+str(i)+",SC"+str(j))
        output.write('\n')
output.write('\n')


output.write('PRINT STRIDE='+str(args_stride)+ ' ARG=* FILE=COLVAR_SC')
output.write('\n')
output.write('\n')
output.close()

print("Running plumed_SC.dat on the folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_SC.dat --mf_xtc {args_folded_trajectory} --pdb {args_reference} --mc {args_mcfile} 1> plumed_SC_folded.out")
print("Running plumed_SC.dat on the unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_SC.dat --mf_xtc {args_unfolded_trajectory} --pdb {args_reference} --mc {args_mcfile} 1> plumed_SC_unfolded.out")
os.chdir(inp_dir)

# read the COLVAR files for the side chain contacts
print("Reading the COLVAR file for the side chains ...")
df_sc = read_colvar('COLVAR_SC', 'contside')
dffilteredFmU = df_sc[df_sc['FminusU'] > 0.0].dropna()
dffilteredFmU = dffilteredFmU[dffilteredFmU['meanF'] > args_cutoff].dropna()
dffilteredFmU = dffilteredFmU[dffilteredFmU['lda'] > args_lda].dropna()
#dffilteredFmU.sort_values(by=['lda'],ascending=False)

dffilteredUmF = df_sc[df_sc['FminusU'] < 0.0].dropna()
dffilteredUmF = dffilteredUmF[dffilteredUmF['meanU'] > args_cutoff].dropna()
dffilteredUmF = dffilteredUmF[dffilteredUmF['lda'] > args_lda].dropna()
#dffilteredUmF.sort_values(by=['lda'],ascending=False)

dffiltered = pd.concat([dffilteredFmU, dffilteredUmF])

##################################################################################
# Adding the solvation bias

# this PLUMED script will compute the solvation for each oxygen, nitrogen and carbon atoms

print("Writing plumed_solvation.dat for the solvation features ...")
output = open('plumed_solvation.dat', 'w')
output.write(f"MOLINFO MOLTYPE=protein STRUCTURE={args_reference_protein}")
output.write('\n')
output.write('WO: GROUP ATOMS='+str(first_oxygen)+'-'+str(last_oxygen)+':'+str(water_model))
output.write('\n')
output.write('WH: GROUP ATOMS='+str(first_oxygen+1)+'-'+str(last_oxygen+1)+':'+str(water_model)+','+str(first_oxygen+2)+'-'+str(last_oxygen+2)+':'+str(water_model))
output.write('\n')

atom_cno = table[(table['element']== 'O') | (table['element']== 'N') | (table['element']== 'C') ].reset_index()

list_atoms = []

for i in range(0,len(atom_cno)):
    serial = atom_cno['serial'].iloc[i]
    atom_type = atom_cno['element'].iloc[i]
    output.write(atom_type+str(serial)+': COORDINATION GROUPA='+str(serial)+ ' GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=10 D_MAX=1.0} #NLIST NL_CUTOFF=1.5 NL_STRIDE=20') ##DEBUG
    output.write('\n')
    list_atoms.append(atom_type+str(serial))

to_write = 'PRINT ARG=' +str(list(list_atoms))+' STRIDE='+str(args_stride)+ ' FILE=COLVAR_solvation'
to_write_clean = to_write.replace("]"," ").replace("["," ").replace(", ",",").replace(" ,",",").replace("= ","=").replace("'","")
output.write(to_write_clean)
output.close()

print("Running plumed_solvation.dat on the folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_solvation.dat --mf_xtc {args_folded_trajectory} --pdb {args_reference} 1> plumed_solvation_folded.out")
print("Running plumed_solvation.dat on the unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_solvation.dat --mf_xtc {args_unfolded_trajectory} --pdb {args_reference} 1> plumed_solvation_unfolded.out")
os.chdir(inp_dir)

os.system(f'sed -i "s/#NLIST/NLIST/g" {inp_dir}/plumed_solvation.dat') ##DEBUG

# read the COLVAR files for the solvation
print("Reading the COLVAR file for the solvation bias ...")
df_solvation = read_colvar('COLVAR_solvation', '')
df_solvation_C = df_solvation[df_solvation['labels'].str.startswith('C')]
# keeping top 4 carbon atoms
df_solvation_C = df_solvation_C.sort_values(by=['lda'],ascending=False).head(4).reset_index()
df_solvation_NO = df_solvation[df_solvation['labels'].str.startswith(('N','O'))]
# keeping top 3 carbon atoms
df_solvation_NO = df_solvation_NO.sort_values(by=['lda'],ascending=False).head(3).reset_index()

############################################################################################
# Writing the final PLUMED file

print("Writing the final plumed file ...")
output = open('plumed_final.dat', 'w')
output.write('# This PLUMED file has been generated by bioinspired_features_generation.py on '+ str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'.')
output.write('\n')
output.write('# Command : python bioinspired_features_generation_v.1.py -F '+ str(args_path_folded) + ' -U ' + str(args_path_unfolded) + ' -r ' + str(args_reference) + ' -rp ' + str(args_reference_protein) + ' -rca ' + str(args_reference_ca) + ' -mc ' + str(args_mcfile) + ' -l ' + str(args_lda) + ' -c ' + str(args_cutoff) + ' -s ' + str(args_stride) + (' -e ' if args_explicit else '') + (' -py ' if args_pymol else '')  + (' -y ' if args_yes else ''))
output.write('\n')
output.write('\n')
output.write("#RESTART")
output.write('\n')
output.write(f"MOLINFO MOLTYPE=protein STRUCTURE={args_reference_protein}")
output.write('\n')
output.write('WO: GROUP ATOMS='+str(first_oxygen)+'-'+str(last_oxygen)+':'+str(water_model))
output.write('\n')
output.write('WH: GROUP ATOMS='+str(first_oxygen+1)+'-'+str(last_oxygen+1)+':'+str(water_model)+','+str(first_oxygen+2)+'-'+str(last_oxygen+2)+':'+str(water_model))
output.write('\n')
output.write(f'rmsd_ca: RMSD REFERENCE={args_reference_ca} TYPE=OPTIMAL') ##TO DEBUG
output.write('\n')
output.write('\n')

#############################################################
# H-bonds ###################################################
#############################################################

for i in range(0,len(df_hard_soft)):
    labels = df_hard_soft['labels'][i]
    match_atomtype = re.findall(r'([A-Za-z0-9]+)_\d+', labels)
    match_number = re.findall(r'_(\d+)', labels)
    original_label = labels.replace("hbond_", "")

    A = match_atomtype[1] # the atom type of A
    A_serial = match_number[1] # the atom number of A
    
    H = match_atomtype[0] # the atom type of H
    H_serial = match_number[0] # the atom number of H
    
    D = hatom_hdonor_dict[match_atomtype[0]] # the atom type of D, corresponding to the match in the dictionnary
    D_resSeq_temp = table[table['serial'] == int(match_number[0])]['resSeq'].reset_index() # temporary variable, dataframe of the atom number of H
    D_resSeq = D_resSeq_temp['resSeq'].loc[0] # residue number of D, which is the same one as H
    D_serial_temp = table[(table['resSeq'] == int(D_resSeq)) & (table['name'] == D)]['serial'].reset_index() # temporary variable, dataframe of the residue corresponding to the residue number of D and its atom type
    D_serial = D_serial_temp['serial'].loc[0] # atom number of D
    
    B = atom_hacceptor_dict[match_atomtype[1]]

    if isinstance(B,str): # meaning it is a string so B is only one value

        B_resSeq_temp = table[table['serial'] == int(match_number[1])]['resSeq'].reset_index()
        B_resSeq = B_resSeq_temp['resSeq'].loc[0]
        B_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B)]['serial'].reset_index()
        B_serial = B_serial_temp['serial'].loc[0]
    
        label_cont = "cont_"+ original_label
        label_HD = "H-D_"+str(H)+"_"+str(H_serial)+"-"+str(D)+"_"+str(D_serial)
        label_HA = "H-A_"+str(H)+"_"+str(H_serial)+"-"+str(A)+"_"+str(A_serial)
        label_HB = "H-B_"+str(H)+"_"+str(H_serial)+"-"+str(B)+"_"+str(B_serial)
        label_DA = "D-A_"+str(D)+"_"+str(D_serial)+"-"+str(A)+"_"+str(A_serial)
        label_AB = "A-B_"+str(A)+"_"+str(A_serial)+"-"+str(B)+"_"+str(B_serial)
        label_ang_DHA = "ang_DHA_"+str(labels)
        label_ang_BAH = "ang_BAH_"+str(labels)
        label_hbond = "hbond_"+str(original_label)
        label_VA = "VA_"+str(original_label)
        label_VH = "VH_"+str(original_label)
        label_NWH = "NWH_"+str(original_label)
        label_NWO = "NWO_"+str(original_label)
    
        s_cont=label_cont+": COORDINATION GROUPA="+str(match_number[0])+" GROUPB="+str(match_number[1])+" SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=8}"
        sHD=label_HD+": DISTANCE ATOMS="+str(H_serial)+","+str(D_serial)
        sHA=label_HA+": DISTANCE ATOMS="+str(H_serial)+","+str(A_serial)
        sHB=label_HB+": DISTANCE ATOMS="+str(H_serial)+","+str(B_serial)
        sDA=label_DA+": DISTANCE ATOMS="+str(D_serial)+","+str(A_serial)
        sAB=label_AB+": DISTANCE ATOMS="+str(A_serial)+","+str(B_serial)
        sang_DHA = label_ang_DHA+": CUSTOM ARG="+label_DA+","+label_HD+","+label_HA+" FUNC=x/(y+z) PERIODIC=NO"
        sang_BAH = label_ang_BAH+": CUSTOM ARG="+label_HB+","+label_HA+","+label_AB+" FUNC=x/(y+z) PERIODIC=NO"
        shbond=label_hbond+": CUSTOM ARG="+label_cont+","+label_ang_DHA+","+label_ang_BAH+" FUNC=x*y*z PERIODIC=NO"
        s_vA = label_VA +": GHOST ATOMS="+str(B_serial)+","+str(A_serial)+","+str(int(A_serial)+1)+" COORDINATES=0.25,0.0,0.0"
        s_NWH = label_NWH +": COORDINATION GROUPA="+ label_VA +" GROUPB=WH SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20"
        s_VH = label_VH +": GHOST ATOMS="+str(D_serial)+","+str(H_serial)+","+str(int(H_serial)+1)+" COORDINATES=0.25,0.0,0.0"
        s_NWO = label_NWO +": COORDINATION GROUPA="+ label_VH +" GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20"
        
        output.write("# Distance between pairs of atoms of residues "+str(D_resSeq)+" and "+str(B_resSeq))
        output.write('\n')
        output.write(s_cont)
        output.write('\n')
        output.write(sHD)
        output.write('\n')
        output.write(sHA)
        output.write('\n')
        output.write(sHB)
        output.write('\n')
        output.write(sDA)
        output.write('\n')
        output.write(sAB)
        output.write('\n')
        output.write(sang_DHA)
        output.write('\n')
        output.write(sang_BAH)
        output.write('\n')
        output.write(shbond)
        output.write('\n')
        output.write(s_vA)
        output.write('\n')
        output.write(s_VH)
        output.write('\n')

        if args_explicit:
            output.write(s_NWH)
            output.write('\n')
            output.write(s_NWO)
            output.write('\n')

    else: # meaning it is ND1 from histidine or SD from methionine that are bound to two different B => COM between two B
        B_resSeq_temp = table[table['serial'] == int(match_number[1])]['resSeq'].reset_index()
        B_resSeq = B_resSeq_temp['resSeq'].loc[0]

        B_1 = B[0]
        B_2 = B[1]
        
        B_1_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B_1)]['serial'].reset_index()
        B_1_serial = B_1_serial_temp['serial'].loc[0]

        B_2_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B_2)]['serial'].reset_index()
        B_2_serial = B_2_serial_temp['serial'].loc[0]


        label_cont = "cont_"+str(df3['labels'][i])
        label_com = "com_"+str(B_1)+"_"+str(B_1_serial)+"-"+str(B_2)+"_"+str(B_2_serial)
        label_HD = "H-D_"+str(H)+"_"+str(H_serial)+"-"+str(D)+"_"+str(D_serial)
        label_HA = "H-A_"+str(H)+"_"+str(H_serial)+"-"+str(A)+"_"+str(A_serial)
        label_HB = "H-B_"+str(H)+"_"+str(H_serial)+"-"+str(B_1)+"_"+str(B_1_serial) # for the sake of the nedt steps I kept only B_1 for the name
        label_DA = "D-A_"+str(D)+"_"+str(D_serial)+"-"+str(A)+"_"+str(A_serial)
        label_AB = "A-B_"+str(A)+"_"+str(A_serial)+"-"+str(B_1)+"_"+str(B_1_serial)
        label_ang_DHA = "ang_DHA_"+str(labels)
        label_ang_BAH = "ang_BAH_"+str(labels)
        label_hbond = "hbond_"+str(labels)
        label_VA = "VA_"+str(original_label)
        label_VH = "VH_"+str(original_label)
        label_NWH = "NWH_"+str(original_label)
        label_NWO = "NWO_"+str(original_label)

        s_cont=label_cont+": COORDINATION GROUPA="+str(match_number[0])+" GROUPB="+str(match_number[1])+" SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=8}"
        com = label_com+": COM ATOMS="+str(B_1_serial)+","+str(B_2_serial)
        sHD=label_HD+": DISTANCE ATOMS="+str(H_serial)+","+str(D_serial)
        sHA=label_HA+": DISTANCE ATOMS="+str(H_serial)+","+str(A_serial)
        sHB=label_HB+": DISTANCE ATOMS="+str(H_serial)+","+str(label_com)
        sDA=label_DA+": DISTANCE ATOMS="+str(D_serial)+","+str(A_serial)
        sAB=label_AB+": DISTANCE ATOMS="+str(A_serial)+","+str(label_com)
        sang_DHA = label_ang_DHA+": CUSTOM ARG="+label_DA+","+label_HD+","+label_HA+" FUNC=x/(y+z) PERIODIC=NO"
        sang_BAH = label_ang_BAH+": CUSTOM ARG="+label_HB+","+label_HA+","+label_AB+" FUNC=x/(y+z) PERIODIC=NO"
        shbond=label_hbond+": CUSTOM ARG="+label_cont+","+label_ang_DHA+","+label_ang_BAH+" FUNC=x*y*z PERIODIC=NO"
        s_VA = label_VA +": GHOST ATOMS="+str(com)+","+str(A_serial)+","+str(int(A_serial)+1)+" COORDINATES=0.25,0.0,0.0"  
        s_VH = label_VH +": GHOST ATOMS="+str(D_serial)+","+str(H_serial)+","+str(int(H_serial)+1)+" COORDINATES=0.25,0.0,0.0"
        s_NWH = label_NWH +": COORDINATION GROUPA="+ label_VA +" GROUPB=WH SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20"
        s_NWO = label_NWO +": COORDINATION GROUPA="+ label_VH +" GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20"
        
        output.write("# Distance between pairs of atoms of residues "+str(D_resSeq)+" and "+str(B_resSeq))
        output.write('\n')
        output.write(s_cont)
        output.write('\n')
        output.write(sHD)
        output.write('\n')
        output.write(sHA)
        output.write('\n')
        output.write(sHB)
        output.write('\n')
        output.write(sDA)
        output.write('\n')
        output.write(sAB)
        output.write('\n')
        output.write(sang_DHA)
        output.write('\n')
        output.write(sang_BAH)
        output.write('\n')
        output.write(shbond)
        output.write('\n')
        output.write(s_vA)
        output.write('\n')
        output.write(s_VH)
        output.write('\n')

        if args_explicit:
            output.write(s_NWH)
            output.write('\n')
            output.write(s_NWO)
            output.write('\n')

############################################

# add the exclusion list based on the ghost atoms

df_VH = df_virtual[df_virtual['labels'].str.contains("virtual_HD")].reset_index()
df_VA = df_virtual[df_virtual['labels'].str.contains("virtual_AB")].reset_index()

HD_list = []
NPA = ""

for i in range(0,len(df_VH)):
    label = df_VH['labels'][i]
    match_atomtype = re.findall(r'([A-Za-z0-9]+)_\d+', label)
    match_number = re.findall(r'_(\d+)', label)
    label_temp = match_atomtype[0] + "_" + str(match_number[0]) + "-" + match_atomtype[1] + "_" + str(match_number[1])
    temp = df_VH[df_VH['labels'].str.contains(label_temp)].reset_index()

    serial = ""
    for j in range(0,len(temp)):
        label = temp['labels'][j]
        match_number = re.findall(r'_(\d+)', label)
        serial +=str(match_number[2])+","
        
    label_VH_coordination = "NPA_" + label_temp
    s = label_VH_coordination + ": COORDINATION GROUPA=" + "VH_" +label_temp + " GROUPB=" + serial +" SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20"

    if s not in HD_list:
       output.write(s)
       output.write('\n')
       HD_list.append(s)
       NPA+=str(label_VH_coordination)+","
    else :
       continue

output.write('\n')

AB_list = []

NPD = ""

for i in range(0,len(df_VA)):
    label = df_VA['labels'][i]
    match_atomtype = re.findall(r'([A-Za-z0-9]+)_\d+', label)
    match_number = re.findall(r'_(\d+)', label)
    label_temp = match_atomtype[0] + "_" + str(match_number[0]) + "-" + match_atomtype[1] + "_" + str(match_number[1])
    temp = df_VA[df_VA['labels'].str.contains(label_temp)].reset_index()

    serial = ""
    for j in range(0,len(temp)):
        label = temp['labels'][j]
        match_number = re.findall(r'_(\d+)', label)
        serial +=str(match_number[2])+","
        
    label_VA_coordination = "NPD_" + label_temp
    
    s = label_VA_coordination + ": COORDINATION GROUPA=" + "VA_" + label_temp + " GROUPB=" + serial +" SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20"

    if s not in AB_list:
        output.write(s)
        output.write('\n')
        NPD+=str(label_VA_coordination)+","
        AB_list.append(s)
    else :
        continue

#####################################################################################################################
    
output.write('\n')
output.write('\n')
output.write('# Speeding up H-bonds')
output.write('\n')

listcv = ""

if df_hard_sigF.empty == False:
    cont_HB_hardF =""
    NWH_hardF =""
    NWO_hardF =""
    NPA_hardF =""
    NPD_hardF =""

    listcv +="cont_HB_hardF,NWH_hardF,NWO_hardF,NPA_hardF,NPD_hardF,"
    
    for i in range(0,len(df_hard_sigF)):
            desc=str(df_hard_sigF['labels'].iloc[i])
            NWH=desc.replace('hbond','VA' )
            NWO=desc.replace('hbond','VH')
            NPA=desc.replace('hbond','NPA')
            NPD=desc.replace('hbond','NPD')
    
            cont_HB_hardF+=str(desc)+","
            NWH_hardF+=str(NWH)+","
            NWO_hardF+=str(NWO)+","
            NPA_hardF+=str(NPA)+","
            NPD_hardF+=str(NPD)+","
    
    output.write('cont_HB_hardF: COMBINE ARG='+cont_HB_hardF+' PERIODIC=NO')
    output.write('\n')
    output.write('NWH_hardF: COORDINATION GROUPA='+NWH_hardF+' GROUPB=WH SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20')
    output.write('\n')
    output.write('NWO_hardF: COORDINATION GROUPA='+NWO_hardF+' GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20')
    output.write('\n')
    output.write('NPA_hardF: COMBINE ARG=' + NPA_hardF + ' PERIODIC=NO')
    output.write('\n')
    output.write('NPD_hardF: COMBINE ARG=' + NPD_hardF + ' PERIODIC=NO')
    output.write('\n')

if df_hard_sigU.empty == False:
    cont_HB_hardU =""
    NWH_hardU =""
    NWO_hardU =""
    NPA_hardU =""
    NPD_hardU =""

    listcv +="cont_HB_hardU,NWH_hardU,NWO_hardU,NPA_hardU,NPD_hardU,"
    
    for i in range(0,len(df_hard_sigU)):
            desc=str(df_hard_sigF['labels'].iloc[i])
            NWH=desc.replace('hbond','VA' )
            NWO=desc.replace('hbond','VH')
            NPA=desc.replace('hbond','NPA')
            NPD=desc.replace('hbond','NPD')
    
            cont_HB_hardU+=str(desc)+","
            NWH_hardU+=str(NWH)+","
            NWO_hardU+=str(NWO)+","
            NPA_hardU+=str(NPA)+","
            NPD_hardU+=str(NPD)+","

    
    output.write('cont_HB_hardU: COMBINE ARG='+cont_HB_hardU+' PERIODIC=NO')
    output.write('\n')
    output.write('NWH_hardU: COORDINATION GROUPA='+NWH_hardU+' GROUPB=WH SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20')
    output.write('\n')
    output.write('NWO_hardU: COORDINATION GROUPA='+NWO_hardU+' GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20')
    output.write('\n')
    output.write('NPA_hardU: COMBINE ARG=' + NPA_hardU + ' PERIODIC=NO')
    output.write('\n')
    output.write('NPD_hardU: COMBINE ARG=' + NPD_hardU + ' PERIODIC=NO') 
    output.write('\n')

if df_soft_sigF.empty == False:
    cont_HB_softFa = ""
    cont_HB_softFb = ""
    NPA_softF =""
    NPD_softF =""

    listcv +="cont_HB_softF,NPA_softF,NPD_softF,"
    
    for i in range(0,len(df_soft_sigF)):
            desc=str(df_soft_sigF['labels'].iloc[i])
            match_number = re.findall(r'_(\d+)', desc)
            NPA=desc.replace('hbond','NPA')
            NPD=desc.replace('hbond','NPD')
    
            ind1 = match_number[0]
            ind2 = match_number[1]
    
            cont_HB_softFa +=str(ind1)+","
            cont_HB_softFb+=str(ind2)+","
            NPA_softF+=str(NPA)+","
            NPD_softF+=str(NPD)+","

    output.write('\n')
    output.write('cont_HB_softF: COORDINATION GROUPA='+cont_HB_softFa+' GROUPB='+cont_HB_softFb +' SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=8} PAIR') 
    output.write('\n')
    output.write('NPA_softF: COMBINE ARG=' + NPA_softF + ' PERIODIC=NO')
    output.write('\n')
    output.write('NPD_softF: COMBINE ARG=' + NPD_softF + ' PERIODIC=NO') 
    output.write('\n')

if df_soft_sigU.empty == False:
    cont_HB_softUa = ""
    cont_HB_softUb = ""
    NPA_softU =""
    NPD_softU =""

    listcv +="cont_HB_softU,NPA_softU,NPD_softU,"
    
    for i in range(0,len(df_soft_sigU)):
            desc=str(df_soft_sigU['labels'].iloc[i])
            match_number = re.findall(r'_(\d+)', desc)
            NPA=desc.replace('hbond','NPA')
            NPD=desc.replace('hbond','NPD')
    
            ind1 = match_number[0]
            ind2 = match_number[1]
    
            cont_HB_softUa +=str(ind1)+","
            cont_HB_softUb+=str(ind2)+","
            NPA_softU+=str(NPA)+","
            NPD_softU+=str(NPD)+","

    output.write('cont_HB_softU: COORDINATION GROUPA='+cont_HB_softUa+' GROUPB='+cont_HB_softUb +' SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=8} PAIR') 
    output.write('\n') 
    output.write('NPA_softU: COMBINE ARG=' + NPA_softU + ' PERIODIC=NO')
    output.write('\n')
    output.write('NPD_softU: COMBINE ARG=' + NPD_softU + ' PERIODIC=NO') 
    output.write('\n')

output.write('\n')
output.write('# Total')
output.write('\n')  
  
coefficients= {
    "cont_HB_hardF": 1.0,
    "NWH_hardF": -1.0/16.0,
    "NWO_hardF": -1.0/8.0,
    "NPA_hardF": -1.0, ##DEBUG
    "NPD_hardF": -0.5, ##DEBUG
    
    "cont_HB_hardU": -1.0,
    "NWH_hardU": 1.0/16.0,
    "NWO_hardU": 1.0/8.0,
    "NPA_hardU": 1.0, ##DEBUG
    "NPD_hardU": 0.5, ##DEBUG
    
    "cont_HB_softF": 1.0,
    "NPA_softF": -1.0, ##DEBUG
    "NPD_softF": -0.5, ##DEBUG
    
    "cont_HB_softU": -1.0,
    "NPA_softU": 0.0,
    "NPD_softU": 0.0,
}

temp = listcv.split(",")
coeff = ""

for i in range(0, len(temp)-1):
    temp2 = coefficients[temp[i]]
    coeff +=str(temp2)+","

label_cmap_H = 'diffHB_compact'
output.write(label_cmap_H +': COMBINE ARG=' + listcv + ' COEFFICIENTS='+ coeff +' PERIODIC=NO')
output.write('\n') 

###########################################
# normal explicit definitions
################################

if args_explicit:
     
    diff_hbond = []

    output.write('\n')
    output.write('# Combine the H-bonds, non-compact')

    for i in range(0,len(df_hard_soft)):
        labels = df_hard_soft['labels'][i]
        original_label = labels.replace("hbond_", "")
        label_final = "diff_" + str(original_label)
        label_NPD = "NPD_"+str(original_label)
        label_NPA = "NPA_"+str(original_label)
        label_NWH = "NWH_"+str(original_label)
        label_NWO = "NWO_"+str(original_label)
        label_hbond = "cont_"+str(original_label)
        label_final = "diff_" + str(original_label) 

        if labels in list(df_hard_sigF['labels']): # for hard folded
            s_label_final = label_final + ": CUSTOM ARG=" + labels  + "," + label_NWH + "," + label_NWO + "," + label_NPA  + "," + label_NPD + " VAR=x,y,z,t,w FUNC=x-y/16-z/8-t-w/2 PERIODIC=NO" ##DEBUG
        if labels in list(df_hard_sigU['labels']): # for hard unfolded
            s_label_final = label_final + ": CUSTOM ARG=" + labels  + "," + label_NWH + "," + label_NWO + "," + label_NPA  + "," + label_NPD + " VAR=x,y,z,t,w FUNC=-x+y/16+z/8+t+w/2 PERIODIC=NO" ##DEBUG
        if labels in list(df_soft_sigF['labels']): # for soft folded
            s_label_final = label_final + ": CUSTOM ARG=" + label_hbond + "," + label_NPA  + "," + label_NPD + " VAR=x,t,w FUNC=x-t-w/2 PERIODIC=NO" ##DEBUG
        if labels in list(df_soft_sigU['labels']): # for soft unfolded
            s_label_final = label_final + ": CUSTOM ARG=" + label_hbond + ", VAR=x FUNC=-x PERIODIC=NO"

        diff_hbond.append(label_final)
        output.write('\n')
        output.write(s_label_final)
        output.write('\n')

    output.write('\n')
    s = "diffHB_non_compact: COMBINE ARG=" + str(list(diff_hbond)) + " COEFFICIENTS=" + str(len(diff_hbond)*"1.0,") + "  PERIODIC=NO"
    s2 = s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","=").replace("'","")
    output.write(s2)

#############################################################
# Side chains ###############################################
#############################################################

output.write('\n')
output.write('\n')
output.write("# Now for the side chains ")
output.write('\n') 

target = mdtraj.load(args_reference_protein)

#to save the amount of residues forming the protein
numres=target.topology.select("name CA").shape

# to calculate the geometric center of the side chains (SC)
resID_nogly=[]
for i in range(0,numres[0]):
    p=target.topology.select("resid "+str(i)+" and not backbone and not (type H)")+1
    p = p.tolist()
    if len(p) !=0:
        s="SC"+str(i+1)+": CENTER ATOMS="+str(p)+"MASS"
        output.write(s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","="))
        output.write('\n')
    else:
        p=target.topology.select("resid "+str(i)+" and backbone and name CA")+1
        s="SC"+str(i+1)+": CENTER ATOMS="+str(p)+"MASS"
        output.write(s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","="))
        output.write('\n')
output.write('\n')

#CV of the distances between the SCs (adjacent SC are not taken into account)
#the list resID_nogly is needed to take into account the removal of GLYs

if args_explicit: 
    for i in range(1,int(numres[0])+1):
        for j in range(i+2,int(numres[0])+1):
            label = "contside"+str(i)+"-"+str(j)+": "
            output.write(str(label)+"COORDINATION GROUPA=SC"+str(i)+" GROUPB=SC"+str(j)+"  SWITCH={RATIONAL D_0=0.0 R_0=0.80 NN=4 MM=8}")
            output.write('\n')
    output.write('\n')

#calculate the point where to evaluate water
#SC2-19: CENTER ATOMS=SC2,SC19
for i in range(1,int(numres[0])+1):
    for j in range(i+2,int(numres[0])+1):
        label = "SC"+str(i)+"-"+str(j)+": "
        output.write(str(label)+"CENTER ATOMS=SC"+str(i)+",SC"+str(j))
        output.write('\n')
output.write('\n')

list_diff = []

listsc = ""

if dffilteredFmU.empty == False:
    contsideaF =""
    contsidebF =""
    exclusion_listF =""

    listsc +="contsideF,exSCF"
    
    for i in range(0,len(dffilteredFmU)):
            desc=str(dffilteredFmU['labels'].iloc[i])
            exclabel=desc.replace("contside", "exc_SC" )
            difflabel=desc.replace("contside", "diff_SC" )
            exclusion_listF+= str(exclabel)+","
            excgroupa=desc.replace("contside", "SC" )
            excindeces=desc.replace("contside", "" )
            ind1 = int(excindeces.split("-")[0])
            ind2 = int(excindeces.split("-")[1])
    
            contsideaF+="SC"+str(ind1)+","
            contsidebF+="SC"+str(ind2)+","
            
            excgroupb=""
            for j in range(1,int(numres[0])+1):
                if (j != ind1) and (j != ind2):
                    excgroupb+="SC"+str(j)+","
            output.write('\n')    
            output.write(exclabel+": COORDINATION GROUPA="+excgroupa+" GROUPB="+excgroupb+" SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20") 
            output.write('\n')

            if args_explicit: 
                output.write(difflabel+": COMBINE ARG="+desc+","+exclabel+" COEFFICIENTS=1.0,-1.0 PERIODIC=NO")
                list_diff.append(difflabel)
                output.write('\n')
        
    output.write('\n')
    output.write("contsideF: COORDINATION GROUPA="+contsideaF+" GROUPB="+contsidebF+" SWITCH={RATIONAL D_0=0.0 R_0=0.80 NN=4 MM=8} PAIR") 
    output.write('\n')
    output.write("exSCF: COMBINE ARG="+exclusion_listF+" PERIODIC=NO")
    output.write('\n')


if dffilteredUmF.empty == False:
    contsideaU =""
    contsidebU =""
    exclusion_listU =""
    
    if dffilteredFmU.empty == False:
        listsc +=",contsideU,exSCU"
    else:
        listsc +="contsideU,exSCU"
    
    for i in range(0,len(dffilteredUmF)):
            desc=str(dffilteredUmF['labels'].iloc[i])
            exclabel=desc.replace("contside", "exc_SC" )
            difflabel=desc.replace("contside", "diff_SC" )
            exclusion_listU+= str(exclabel)+","
            excgroupa=desc.replace("contside", "SC" )
            excindeces=desc.replace("contside", "" )
            ind1 = int(excindeces.split("-")[0])
            ind2 = int(excindeces.split("-")[1])
    
            contsideaU+="SC"+str(ind1)+","
            contsidebU+="SC"+str(ind2)+","
            
            excgroupb=""
            for j in range(1,int(numres[0])+1):
                if (j != ind1) and (j != ind2):
                    excgroupb+="SC"+str(j)+","
            output.write('\n')    
            output.write(exclabel+": COORDINATION GROUPA="+excgroupa+" GROUPB="+excgroupb+" SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20") 
            output.write('\n')

            if args_explicit:     
                output.write(difflabel+": COMBINE ARG="+desc+","+exclabel+" COEFFICIENTS=-1.0,1.0 PERIODIC=NO")
                output.write('\n')
                list_diff.append(difflabel)


    output.write('\n')
    output.write("contsideU: COORDINATION GROUPA="+contsideaU+" GROUPB="+contsidebU+" SWITCH={RATIONAL D_0=0.0 R_0=0.80 NN=4 MM=8} PAIR") 
    output.write('\n')
    output.write("exSCU: COMBINE ARG="+exclusion_listU+" PERIODIC=NO") 
    output.write('\n')


temp = listsc.split(",")

coefficients_SC= {
    "contsideF": 1.0,
    "exSCF": -1.0,
    "contsideU": -1.0,
    "exSCU": 1.0
}

coeff_SC = ""

for i in range(0, len(temp)):
    temp2 = coefficients_SC[temp[i]]
    coeff_SC +=str(temp2)+","

label_cmap_SC = 'cmap_compact'
output.write(label_cmap_SC+': COMBINE ARG=' + listsc + ' COEFFICIENTS='+ coeff_SC +' PERIODIC=NO')
output.write('\n')

if args_explicit:
    output.write('\n')
    output.write('# Combine the side chain cmap')
    output.write('\n')
    output.write('\n')
    label_cmap_SC_non_compact = "cmap_non_compact"
    s_cmap = label_cmap_SC_non_compact + ": COMBINE ARG=" + str(list(list_diff)) + " COEFFICIENTS=" + str(len(list_diff)*"1.0,") + "  PERIODIC=NO"
    s2_cmap = s_cmap.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","=").replace("'","")
    output.write(s2_cmap)

###############################################################################
# Adding the coordination to water as auxiliary CV

# add finally the carbon/oxygen/nitrogen atoms coordination to water

list_atoms = []

for i in range(0,len(df_solvation_C)):
    label = df_solvation_C['labels'][i]
    serial = re.findall(r'(\d+)', label)

    output.write('\n')
    output.write(str(label)+': COORDINATION GROUPA='+str(serial[0])+ ' GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=10 D_MAX=1.0} NLIST NL_CUTOFF=1.5 NL_STRIDE=20')
    output.write('\n')
    list_atoms.append(str(label))
    
for i in range(0,len(df_solvation_NO)):
    label = df_solvation_NO['labels'][i]
    serial = re.findall(r'(\d+)', label)
    
    output.write(str(label)+': COORDINATION GROUPA='+str(serial[0])+ ' GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=10 D_MAX=1.0} NLIST NL_CUTOFF=1.5 NL_STRIDE=20')
    output.write('\n')
    list_atoms.append(str(label))

output.write('\n')

if args_explicit:
    to_write = 'PRINT ARG=diffHB_compact,cmap_compact,diffHB_non_compact,cmap_non_compact STRIDE='+str(stride)+ ' FILE=COLVAR_diff'
else:
        to_write = 'PRINT ARG=diffHB_compact,cmap_compact STRIDE='+str(args_stride)+ ' FILE=COLVAR_diff'

to_write_clean = to_write.replace("]"," ").replace("["," ").replace(", ",",").replace(" ,",",").replace("= ","=").replace("'","")
output.write(to_write_clean)
output.close()

# to remove duplicate rows that would mess up with PLUMED
lines_seen = set() # holds lines already seen
outfile = open('plumed_final_noduplicate.dat', "w")
for line in open('plumed_final.dat', "r"):
    if line not in lines_seen: # not a duplicate
        outfile.write(line)
        lines_seen.add(line)
outfile.close()

print("Running the final plumed file on the folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f'sed -i "s/NLIST/#NLIST/g" {inp_dir}/plumed_final_noduplicate.dat') ##DEBUG
os.system(f"plumed driver --plumed {inp_dir}/plumed_final_noduplicate.dat --mf_xtc {args_folded_trajectory} --pdb {args_reference} --mc {args_mcfile} 1> plumed_final_folded.out")
print("Running the final plumed file on the unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_final_noduplicate.dat --mf_xtc {args_unfolded_trajectory} --pdb {args_reference} --mc {args_mcfile} 1> plumed_final_unfolded.out")
os.chdir(inp_dir)
os.system(f'sed -i "s/#NLIST/NLIST/g" {inp_dir}/plumed_final_noduplicate.dat') ##DEBUG

print("Reading the final COLVAR file ...")
df_final = read_colvar('COLVAR_diff', '')

nb_frames = len(pd.read_csv(args_folded_dir+'/COLVAR_diff', sep=r'\s+',skiprows=1, header=None)-1)


# Writing final results to log file
logging.info(' With a stride of ' + str(int(args_stride)) + ', the number of frames used for the filtering was ' + str(nb_frames))
logging.info(' We recommend to include at least 1000 frames to obtain good statistics')
logging.info('\n')
logging.info(' Number of hard H-bonds significant in folded state: '+ str(len(df_hard_sigF)))
logging.info(' Number of hard H-bonds significant in unfolded state: '+ str(len(df_hard_sigU)))
logging.info(' Number of soft H-bonds significant in folded state: '+ str(len(df_soft_sigF)))
logging.info(' Number of soft H-bonds significant in unfolded state: '+ str(len(df_soft_sigU)))
logging.info(' Number of side chain contacts significant in folded state: '+ str(len(dffilteredFmU)))
logging.info(' Number of side chain contacts significant in unfolded state: '+ str(len(dffilteredUmF)))
logging.info('\n')
for label, meanF, meanU in zip(df_final['labels'], df_final['meanF'], df_final['meanU']):
    logging.info(f" meanF for {label}: {round(meanF, 2)}")
    logging.info(f" meanU for {label}: {round(meanU, 2)}")
logging.info(' Finished on '+ str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'.')
logging.info('\n')

# output pymol session if requested
if args_pymol:
    # load the protein and show it as licorice
    cmd.load(args_reference_protein)
    cmd.show_as("licorice", "all")

    # loop on Hard-Folded, Hard-Unfolded, Soft-Folded, Soft-Unfolded hbonds
    for index, row in df_hard_soft.iterrows():
        # get atoms
        label_name = row['labels'].split("_")
        atom_a = label_name[2].split("-")[0]
        atom_b = label_name[3]
        # get type of hbond and prepare label name
        strength = "Hard" if row['Hard'] else "Soft"
        state = "Folded" if row ['Significant folded'] else "Unfolded"
        dist_name = f"hbond_{strength}_{state}" + "_lda_{0:.2f}".format(row['lda'])
        # draw the distance, hide the label (confusing), and color accordingly
        cmd.distance(dist_name, f"id {atom_a}", f"id {atom_b}")
        cmd.hide("labels", dist_name)
        cmd.color("red" if row['Hard'] else "yelloworange", dist_name)
 
    # loop on folded contacts
    for index, row in (dffilteredFmU.sort_values(by=['lda'],ascending=False)).iterrows():
        # get atoms
        label_name = row['labels'].split("-")
        res_a = ''.join(char for char in label_name[0] if char.isdigit())
        res_b = label_name[1]
        cont_name = f"contact_folded_res{res_a}_res{res_b}" + "_lda_{0:.2f}".format(row['lda'])
        # select two residues sidechains and put a pseudo-atom in the centre of mass
        cmd.select("side_a", f"resi {res_a} and sidechain")
        cmd.select("side_b", f"resi {res_b} and sidechain")
        com_a = cmd.centerofmass("side_a")
        com_b = cmd.centerofmass("side_b")
        cmd.pseudoatom("com_a", pos=com_a)
        cmd.pseudoatom("com_b", pos=com_b)
        # draw a blue distance for folded main contacts and get rid of the rest
        cmd.distance(cont_name, "com_a", "com_b")
        cmd.hide("labels", cont_name)
        cmd.color("blue", cont_name)
        cmd.delete("side_a")
        cmd.delete("side_b")
        cmd.delete("com_a")
        cmd.delete("com_b")
        cmd.hide("labels", cont_name)
        cmd.color("blue", cont_name)

    # loop on unfolded contacts
    for index, row in (dffilteredUmF.sort_values(by=['lda'],ascending=False)).iterrows():
        # get atoms
        label_name = row['labels'].split("-")
        res_a = ''.join(char for char in label_name[0] if char.isdigit())
        res_b = label_name[1]
        cont_name = f"contact_unfolded_res{res_a}_res{res_b}" + "_lda_{0:.2f}".format(row['lda'])
        # select two residues sidechains and put a pseudo-atom in the centre of mass
        cmd.select("side_a", f"resi {res_a} and sidechain")
        cmd.select("side_b", f"resi {res_b} and sidechain")
        com_a = cmd.centerofmass("side_a")
        com_b = cmd.centerofmass("side_b")
        cmd.pseudoatom("com_a", pos=com_a)
        cmd.pseudoatom("com_b", pos=com_b)
        # draw a light-blue distance for folded main contacts and get rid of the rest
        cmd.distance(cont_name, "com_a", "com_b")
        cmd.hide("labels", cont_name)
        cmd.color("lightblue", cont_name)
        cmd.delete("side_a")
        cmd.delete("side_b")
        cmd.delete("com_a")
        cmd.delete("com_b")

    # save session and quit
    print(f"Saving hbond pymol session in {inp_dir}/summary_pymol_session.pse")
    cmd.save("summary_pymol_session.pse")

print('Done ! Check the bioinspired_features.log for details.')
