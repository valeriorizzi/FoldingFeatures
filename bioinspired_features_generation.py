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

version = 2.0
error = '--- ERROR: %s \n'
logging.basicConfig(filename='bioinspired_features.log', filemode='a', format='%(levelname)s:%(message)s', level=logging.INFO)

colours = datetime.datetime.now().strftime("%H:%M:%S:%f").split(":")

_p = "\033[38;5;" + str(int((float(colours[0])*float(colours[3])) % 256)) + "m+\033[0m"
_m = "\033[38;5;" + str(int((float(colours[1])*float(colours[3])) % 256)) + "m-\033[0m"
_h = "\033[38;5;" + str(int((float(colours[2])*float(colours[3])) % 256)) + "m#\033[0m"

logo = "\n"
logo += r"______ _       _                 _              _       " + f"           {_p}{_h}{_p}{_m}{_m}{_m}{_p}{_h}                             " + "\n"
logo += r"| ___ (_)     (_)               (_)            | |      " + f"           {_p} {_m}{_m}{_m} {_p}                              " + "\n"
logo += r"| |_/ /_  ___  _ _ __  ___ _ __  _ _ __ ___  __| |      " + f"           {_p} {_p}{_p}                                 " + "\n"
logo += r"| ___ \ |/ _ \| | '_ \/ __| '_ \| | '__/ _ \/ _` |      " + f"         {_p}{_p}  {_p}{_p}       {_p}{_p}                        " + "\n"
logo += r"| |_/ / | (_) | | | | \__ \ |_) | | | |  __/ (_| |      " + f"       {_m}{_p}     {_h}{_h}{_p}  {_m}{_p}   {_h}                       " + "\n"
logo += r"\____/|_|\___/|_|_| |_|___/ .__/|_|_|  \___|\__,_|      " + f"       {_p}               {_p}{_m}{_m}{_m}{_p}                    " + "\n"
logo += r"                          | |                           " + f"        {_p}               {_p}{_m}{_m}{_m}                    " + "\n"
logo += r"                          |_|                           " + f"        {_p}{_p}         {_m}{_m}{_m}{_m}{_m}{_p}{_p}{_p}{_h}                    " + "\n"
logo += r"______         _                                        " + f"         {_p}        {_h}{_p}{_p}{_p}{_p}{_p}{_p}{_h}                      " + "\n"
logo += r"|  ___|       | |                                       " + f"         {_p}          {_p}{_p}{_p}{_m}{_m}                       " + "\n"
logo += r"| |_ ___  __ _| |_ _   _ _ __ ___  ___                  " + f"        {_p}             {_p}{_m}{_m}{_p}                      " + "\n"
logo += r"|  _/ _ \/ _` | __| | | | '__/ _ \/ __|                 " + f"      {_p}{_p}        {_m}{_m}{_m}{_m}{_m}{_p}{_p}{_p}{_h}                       " + "\n"
logo += r"| ||  __/ (_| | |_| |_| | | |  __/\__ \                 " + f"   {_m}{_p}           {_h}{_p}{_p}{_p}{_p}{_p}{_h}{_h}                        " + "\n"
logo += r"\_| \___|\__,_|\__|\__,_|_|  \___||___/                 " + f"                {_h}{_p}{_m}{_m}{_m}{_p}                          " + "\n"
logo += r"                                                        " + f"                  {_p}{_m}{_m}{_m}{_p}                         " + "\n"
logo += r"                                                        " + f"                      {_p}{_p}                        " + "\n"

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
    Returns a dataframe with some statistics on loaded CVs
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
                    
def process_virtual(df_virtual, dfs_sigHB, output, prefix_HD="NPA", prefix_AB="NPD"):
    """
    Takes as input a dataframe df_virtual and list of dataframes df_sigHB
    Computes the NPD and NPA exclusion lists based on the acceptor/donor candidates in df_virtual and
    the significant hard/soft folded/unfolded hbonds in the dfs_sigHB list of dataframes
    Prints them in the new and quicker CoordinationMap framework
    """
    df_VH = df_virtual[df_virtual['labels'].str.contains("virtual_HD")].reset_index()
    df_VA = df_virtual[df_virtual['labels'].str.contains("virtual_AB")].reset_index()

    HD_list = []
    NPA = ""
    k=1 #total counter
    weight = 1.0 #for oxygen in groupB, weight 1

    for i in range(len(df_VH)):
        
        label = df_VH['labels'][i]
        match_atomtype = re.findall(r'([A-Za-z0-9]+)_\d+', label)
        match_number = re.findall(r'_(\d+)', label)
        hbond = match_atomtype[0] + "_" + str(match_number[0]) + "-" + match_atomtype[1] + "_" + str(match_number[1])
        
        # only keep the virtual contacts that are part of one of the significant H-bonds dfs
        # NOTE: each label_temp should be found only in one of the dfs in dfs_sigHB

        for df_sigHB in dfs_sigHB:
            if df_sigHB['labels'].str.contains(hbond).any() == True:
                
                temp = df_VH[df_VH['labels'].str.contains('_' + hbond + '-')].reset_index() # BEWARE OF THE STRING, DANGEROUS FIX

                serial = ""
                for j in range(len(temp)):
                    label = temp['labels'][j]
                    match_number = re.findall(r'_(\d+)', label)
                    serial += str(match_number[2]) + ","

                label_VH_coordination = f"{prefix_HD}_" + hbond
                s2 = (" GROUPA=" + "VH_" + hbond + " GROUPB=" + serial) #simplified s, to create a library
                if k==1:
                    s = (" GROUPA"+str(k)+"=" + "VH_" + hbond + " GROUPB"+str(k)+"=" + serial + " WEIGHT"+str(k)+"="+str(weight) + " SWITCH"+str(k)+"={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} ")
                else:
                    s = (" GROUPA"+str(k)+"=" + "VH_" + hbond + " GROUPB"+str(k)+"=" + serial + " WEIGHT"+str(k)+"="+str(weight))


                if s2 not in HD_list:
                    output.write(s + "\n")
                    k=k+1
                    HD_list.append(s2)
                    NPA += label_VH_coordination + ","

    AB_list = []
    NPD = ""
    weight = 0.5 #for hydrogen in groupB, weight 1

    for i in range(len(df_VA)):
        label = df_VA['labels'][i]
        match_atomtype = re.findall(r'([A-Za-z0-9]+)_\d+', label)
        match_number = re.findall(r'_(\d+)', label)
        hbond = match_atomtype[0] + "_" + str(match_number[0]) + "-" + match_atomtype[1] + "_" + str(match_number[1])
        
        for df_sigHB in dfs_sigHB:
            if df_sigHB['labels'].str.contains(hbond).any() == True:
        
                temp = df_VA[df_VA['labels'].str.contains('_' + hbond + '-')].reset_index() # BEWARE OF THE STRING, DANGEROUS FIX
                serial = ""
                for j in range(len(temp)):
                    label = temp['labels'][j]
                    match_number = re.findall(r'_(\d+)', label)
                    serial += str(match_number[2]) + ","

                label_VA_coordination = f"{prefix_AB}_" + hbond
                s2 = (" GROUPA=" + "VA_" + hbond + " GROUPB=" + serial) #simplified s, to create a library
                s = (" GROUPA"+str(k)+"=" + "VA_" + hbond + " GROUPB"+str(k)+"=" + serial + " WEIGHT"+str(k)+"="+str(weight)) #to be printed

                if s2 not in AB_list:
                    output.write(s + "\n")
                    k=k+1
                    AB_list.append(s2)
                    NPD += label_VA_coordination + ","
                    
#NEW
def process_SC(dfF,dfU,output,number_residues,purpose):
        k=1 #counter of SC contacts
        weight = 1.0 #weight for the SC F contacts

        for i in range(0,len(dfF)):
            desc=str(dfF['labels'].iloc[i])
            excindeces=desc.replace("contside", "" )
            ind1 = int(excindeces.split("-")[0])
            ind2 = int(excindeces.split("-")[1])

            if purpose==1:
                s = "SC"+str(ind1)+"-"+str(ind2)+": CENTER ATOMS="+"SC"+str(ind1)+","+"SC"+str(ind2)
                output.write(s)
                output.write('\n')

            if purpose==2:
                s_cont="   ATOMS"+str(k)+"="+"SC"+str(ind1)+","+"SC"+str(ind2)+" WEIGHT"+str(k)+"="+str(weight)+" SWITCH"+str(k)+"={RATIONAL D_0=0.0 R_0=0.80 NN=4 MM=8}"
                output.write(s_cont)
                output.write('\n')
                k=k+1

            if purpose==3:
                groupb=""
                for j in number_residues:
                    if j==ind1 or j==ind2:
                        donothing=1
                    else:
                        groupb+="SC"+str(j)+","
                if i==0:
                    s_cont="   GROUPA"+str(i+1)+"="+"SC"+str(ind1)+"-"+str(ind2)+" GROUPB"+str(i+1)+"="+str(groupb)+" SWITCH"+str(i+1)+"={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5}"
                else:
                    s_cont="   GROUPA"+str(i+1)+"="+"SC"+str(ind1)+"-"+str(ind2)+" GROUPB"+str(i+1)+"="+str(groupb)
                output.write(s_cont)
                output.write('\n')

        weight = -1.0 #weight for the SC U contacts
        for i in range(0,len(dfU)):
            desc=str(dfU['labels'].iloc[i])
            excindeces=desc.replace("contside", "" )
            ind1 = int(excindeces.split("-")[0])
            ind2 = int(excindeces.split("-")[1])

            if purpose==1:
                s = "SC"+str(ind1)+"-"+str(ind2)+": CENTER ATOMS="+"SC"+str(ind1)+","+"SC"+str(ind2)
                output.write(s)
                output.write('\n')

            if purpose==2:
                s_cont="   ATOMS"+str(k)+"="+"SC"+str(ind1)+","+"SC"+str(ind2)+" WEIGHT"+str(k)+"="+str(weight)+" SWITCH"+str(k)+"={RATIONAL D_0=0.0 R_0=0.80 NN=4 MM=8}"
                output.write(s_cont)
                output.write('\n')
                k=k+1

            if purpose==4:
                groupb=""
                for j in number_residues:
                    if j==ind1 or j==ind2:
                        donothing=1
                    else:
                        groupb+="SC"+str(j)+","
                if i==0:
                    s_cont="   GROUPA"+str(i+1)+"="+"SC"+str(ind1)+"-"+str(ind2)+" GROUPB"+str(i+1)+"="+str(groupb)+" SWITCH"+str(i+1)+"={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5}"
                else:
                    s_cont="   GROUPA"+str(i+1)+"="+"SC"+str(ind1)+"-"+str(ind2)+" GROUPB"+str(i+1)+"="+str(groupb)
                output.write(s_cont)
                output.write('\n')

def process_HB(df_hard_soft, df_hard_sigF,df_hard_sigU,df_soft_sigF,df_soft_sigU, output, purpose):
    df_hard_soft = df_hard_soft.sort_values(by=['lda'],ascending=False)
    s = ""
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

        if isinstance(B,str):
            B_resSeq_temp = table[table['serial'] == int(match_number[1])]['resSeq'].reset_index()
            B_resSeq = B_resSeq_temp['resSeq'].loc[0]
            B_serial_temp = table[(table['resSeq'] == int(B_resSeq)) & (table['name'] == B)]['serial'].reset_index()
            B_serial = B_serial_temp['serial'].loc[0]

            if purpose==1:
                label_VA = "VA_"+str(original_label)
                label_VH = "VH_"+str(original_label)
                s_vA = label_VA +": LINEARGHOST ATOMS="+str(B_serial)+","+str(A_serial)+" COORDINATES=0.25"
                s_VH = label_VH +": LINEARGHOST ATOMS="+str(D_serial)+","+str(H_serial)+" COORDINATES=0.25"
                output.write(s_vA)
                output.write('\n')
                output.write(s_VH)
                output.write('\n')
            elif purpose==2:
                labels = df_hard_soft['labels'][i] #defines the HB as H/F and hard/soft
                if labels in list(df_hard_sigF['labels']): # for hard folded
                    angolo=1
                    weight=1
                if labels in list(df_hard_sigU['labels']): # for hard unfolded
                    angolo=1
                    weight=-1
                if labels in list(df_soft_sigF['labels']): # for soft folded
                    angolo=0
                    weight=1
                if labels in list(df_soft_sigU['labels']): # for soft unfolded
                    angolo=0
                    weight=-1

                if i==0:
                    s_cont="   ATOMS"+str(i+1)+"="+str(H_serial)+","+str(A_serial)+","+str(D_serial)+","+str(B_serial)+" ANGLE"+str(i+1)+"="+str(angolo)+" WEIGHT"+str(i+1)+"="+str(weight)+" SWITCH"+str(i+1)+"={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=8}"
                    output.write(s_cont)
                    output.write('\n')
                else:
                    s_cont="   ATOMS"+str(i+1)+"="+str(H_serial)+","+str(A_serial)+","+str(D_serial)+","+str(B_serial)+" ANGLE"+str(i+1)+"="+str(angolo)+" WEIGHT"+str(i+1)+"="+str(weight)
                    output.write(s_cont)
                    output.write('\n')
            else:
                label_VA = "VA_"+str(original_label)
                label_VH = "VH_"+str(original_label)
                labels = df_hard_soft['labels'][i] #defines the HB as H/F and hard/soft

                if purpose==3: #for F , VH-WO
                    if labels in list(df_hard_sigF['labels']): # for hard folded
                        s += label_VH +","
                elif purpose==4: #for F , VA-WH
                    if labels in list(df_hard_sigF['labels']): # for hard folded
                        s += label_VA +","
                if purpose==5: #for U , VH-WO
                    if labels in list(df_hard_sigU['labels']): # for hard folded
                        s += label_VH +","
                elif purpose==6: #for U , VA-WH
                    if labels in list(df_hard_sigU['labels']): # for hard folded
                        s += label_VA +","

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

            if purpose==1:
                com = label_com+": COM ATOMS="+str(B_1_serial)+","+str(B_2_serial)
                output.write(com)
                output.write('\n')
                label_VA = "VA_"+str(original_label)
                label_VH = "VH_"+str(original_label)
                s_vA = label_VA +": LINEARGHOST ATOMS="+str(label_com)+","+str(A_serial)+" COORDINATES=0.25"
                s_VH = label_VH +": LINEARGHOST ATOMS="+str(D_serial)+","+str(H_serial)+" COORDINATES=0.25"
                output.write(s_vA)
                output.write('\n')
                output.write(s_VH)
                output.write('\n')
            elif purpose==2:
                labels = df_hard_soft['labels'][i] #defines the HB as H/F and hard/soft
                if labels in list(df_hard_sigF['labels']): # for hard folded
                    angolo=1
                    weight=1
                if labels in list(df_hard_sigU['labels']): # for hard unfolded
                    angolo=1
                    weight=-1
                if labels in list(df_soft_sigF['labels']): # for soft folded
                    angolo=0
                    weight=1
                if labels in list(df_soft_sigU['labels']): # for soft unfolded
                    angolo=0
                    weight=-1

                if i==0:
                    s_cont="   ATOMS"+str(i+1)+"="+str(H_serial)+","+str(A_serial)+","+str(D_serial)+","+str(label_com)+" ANGLE"+str(i+1)+"="+str(angolo)+" WEIGHT"+str(i+1)+"="+str(weight)+" SWITCH"+str(i+1)+"={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=8}"
                    output.write(s_cont)
                    output.write('\n')
                else:
                    s_cont="   ATOMS"+str(i+1)+"="+str(H_serial)+","+str(A_serial)+","+str(D_serial)+","+str(label_com)+" ANGLE"+str(i+1)+"="+str(angolo)+" WEIGHT"+str(i+1)+"="+str(weight)
                    output.write(s_cont)
                    output.write('\n')
            else:
                label_VA = "VA_"+str(original_label)
                label_VH = "VH_"+str(original_label)
                labels = df_hard_soft['labels'][i] #defines the HB as H/F and hard/soft

                if purpose==3: #for F , VH-WO
                    if labels in list(df_hard_sigF['labels']): # for hard folded
                        s += label_VH +","
                elif purpose==4: #for F , VA-WH
                    if labels in list(df_hard_sigF['labels']): # for hard folded
                        s += label_VA +","
                if purpose==5: #for U , VH-WO
                    if labels in list(df_hard_sigU['labels']): # for hard folded
                        s += label_VH +","
                elif purpose==6: #for U , VA-WH
                    if labels in list(df_hard_sigU['labels']): # for hard folded
                        s += label_VA +","

    if purpose==3 or purpose==5:
        s += " GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20"
        output.write(s)
    if purpose==4 or purpose==6:
        s += " GROUPB=WH SWITCH={RATIONAL D_0=0.0 R_0=0.35 NN=2 MM=10 D_MAX=0.5} NLIST NL_CUTOFF=0.8 NL_STRIDE=20"
        output.write(s)

# Parse user inputs
parser = argparse.ArgumentParser(
    prog = "bioinspired_features",
    description = "Designs and filters bioinspired features for use in enhanced sampling simulations of peptide folding.",
    epilog = "Thanks for using %(prog)s! ",
)

# compulsory
parser.add_argument("-f", "--folded", required=True, type=os.path.abspath, help="unbiased folded trajectory (xtc file)")
parser.add_argument("-u", "--unfolded", required=True, type=os.path.abspath, help="unbiased unfolded trajectory (xtc file)")
parser.add_argument("-r", "--reference_protein", required=True, type=os.path.abspath, help="reference PDB file for the full system")
parser.add_argument("-mc", "--mcfile", required=True, type=os.path.abspath, help="PLUMED mcfile")

# other options
parser.add_argument("-w", "--where", required=False, default="protein", type=str, help="Selection in mdtraj lingo. Default is 'protein'. If you need a specific part of the protein you can use e.g. 'protein and (residue 6 to 25 or residue 80 to 95)'")
parser.add_argument("-l", "--lda", required=False, default=0.3, type=float, help="lda value to use for the filtering. Default: 0.3")
parser.add_argument("-c", "--cutoff", required=False, default=0.6, type=float, help="cutoff value to use for the filtering. Default: 0.6nm")
parser.add_argument("-s", "--stride", required=False, default=10, type=int, help="STRIDE for the PLUMED file in the filtering steps. Default: 10")
parser.add_argument("-py", "--pymol", required=False, action='store_true', default=False, help='saves a session of pymol with the main hydrogen bonds highlighted; requires the pymol module')
parser.add_argument("-sy", "--symm", required=False, action='store_true', default=False, help='treats the two F and U symmetrically, i.e. it subtracts NNC from the HB_U')

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
args_reference_protein = args.reference_protein
args_mcfile = args.mcfile
# additional parameters
args_where = args.where
args_lda = args.lda
args_cutoff = args.cutoff
args_stride = args.stride
args_pymol = args.pymol
args_symm = args.symm

full_command = 'python bioinspired_features_generation.py -f '+ str(args_path_folded) + ' -u ' + str(args_path_unfolded) + ' -r ' + str(args_reference_protein) \
             + ' -mc ' + str(args_mcfile) + ' -l ' + str(args_lda) + ' -c ' + str(args_cutoff) + ' -s ' + str(args_stride) \
             + (' -py' if args_pymol else '') + (' -sy' if args_symm else '')

# Check format of input
if (args_folded_trajectory[-3:] != "xtc" or args_unfolded_trajectory[-3:] != "xtc"):
    sys.exit(error%('Trajectories must be an xtc file.'))
if args_reference_protein[-3:] != "pdb":
    sys.exit(error%('Reference structures must be a pdb file.'))
# Check files availability
missing_files = []
for infile in [args_reference_protein, args_mcfile, args_path_folded, args_path_unfolded]:
    if not os.path.isfile(infile):
        missing_files.append(infile)
if len(missing_files) > 0:
    sys.exit(error%('The following input files have not been found:\n\n  > ' + '\n  > '.join(missing_files)))
# check plumed sourcing
try:
    subprocess.run(["plumed", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
except:
    sys.exit(error%('PLUMED is not installed or has not been sourced properly.'))

# load reference PDB file (only protein)
full_topology = mdtraj.load(args_reference_protein, top=args_reference_protein).topology
topology = mdtraj.load(args_reference_protein, top=args_reference_protein, atom_indices=full_topology.select(args_where)).topology
number_residues = [res.resSeq for res in topology.residues]
table, bonds = topology.to_dataframe()

# Save the protein locally as pdb file to help for analysis
outname_protein_ref = 'ref_protein_biofeat.pdb'
ref_prot_topology = mdtraj.load(args_reference_protein, top=args_reference_protein, atom_indices=full_topology.select("protein"))
ref_prot_topology[0].save(outname_protein_ref)

# Save the CAs of selection of protein locally as pdb file to help for analysis
outname_protein_CA_ref = 'ref_proteinCA_biofeat.pdb'
ref_prot_CA_topology = mdtraj.load(args_reference_protein, top=args_reference_protein, atom_indices=full_topology.select("name CA and (" + args_where + ")"))
ref_prot_CA_topology[0].save(outname_protein_CA_ref)

# Fix the pdb by cleaning the header and foooter which are not liked by PLUMED
with open(outname_protein_CA_ref, 'r') as f:
    lines_CA = f.readlines()

with open(outname_protein_CA_ref, 'w') as f:
    for line in lines_CA:
        if line.startswith(('ATOM', 'HETATM')):
            f.write(line)

# checks the type of water (# atoms in the model)
topology_solv = mdtraj.load(args_reference_protein, top=args_reference_protein, atom_indices=full_topology.select("water")).topology
ref_solv, bonds_solv = topology_solv.to_dataframe()

n_atoms_water = len(full_topology.select("water"))
n_res_waters = topology_solv.n_residues
water_model = n_atoms_water/n_res_waters

# checks that the number is divisible by three or four, not perfect but relatively robust
if (water_model == 3.0 or water_model == 4.0): # 3 and 4 points supported
    water_model = int(water_model)
elif (water_model % 1 == 0):
    sys.exit(error%('The water model detected has ' + str(water_model) + ' atoms per residue, which is not supported'))
else:
    sys.exit(error%('Non-integer number of water atoms detected (' + str(water_model) + ' atoms per residue). Check the reference file: ' + args_reference_protein))

first_oxygen = list(ref_solv[(ref_solv['resName'] == 'HOH') & (ref_solv['name'] == 'O')].head(1)['serial'])[0]
last_oxygen = first_oxygen + water_model * (n_res_waters - 1)

# check if pymol works if the saving of a pymol session is requested
if args_pymol:
    try:
        import pymol
    except ImportError:
        print("\nWARNING: Module pymol not found. Is it installed?")
        print("WARNING: Will not print the pymol session...")
        args_pymol = False
    except:
        print("WARNING: Cannot use pymol module, but it was found. Something else went wrong!")
        print("WARNING: Will not print the pymol session...")
        args_pymol = False
    else:
        from pymol import cmd

# some on screen info
print("\n")
print(logo)
print("\n")
print("################################################################################################################")
print("Working directory           : " + str(inp_dir))
print("Reference box               : " + str(args_reference_protein))
print("Folded trajectory file      : " + str(args_path_folded))
print("Unfolded trajectory file    : " + str(args_path_unfolded))
print("Reference mcfile            : " + str(args_mcfile))
print("Water model detected        : " + str(water_model) + "-points water model.")
print("Selection for feature ext.  : " + str(args_where))
print("Stride                      : " + str(int(args_stride)))
print("LDA cutoff value            : " + str(args_lda))
print("Cutoff value for the H-bonds: " + str(args_cutoff))
print("Symmetrical F and U         : " + str(args_symm))
print("Pymol session generation    : " + str(args_pymol))
if args_pymol: print(f"Pymol session saved in file : {inp_dir}/symmary_pymol_session.pse")
print("################################################################################################################")
print("\n")

logging.info(' BIOINSPIRED_FEATURES_GENERATION Python script version ' + str(version))
logging.info(' Please read and cite the following reference:')
logging.info(' The Arch from the Stones: Understanding Protein Folding Energy Landscapes via Bioinspired Collective Variables')
logging.info(' The Journal of Physical Chemistry Letters, 16, 9636-9645, (2025)')
logging.info(' Rizzi, V., Héritier, M., Piasentin, N., Aureli, S., & Gervasio, F. L.')
logging.info(' doi: 10.1021/acs.jpclett.5c02079')
logging.info(' Command Line:')
logging.info(' ' + full_command)
logging.info(' Started on '+ str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'.')
logging.info('\n')

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
output.write('# all distances between possible H-bonds')
output.write('\n')

# print out all the combination of possible intra-protein hbonds
for i in range(0, len(list_hbond_donors)): # to go through the H-bond donors list
    for j in range(0, list_hbond_donors[i].shape[0]): 
        for z in range(0, len(list_hbond_acceptors)):
            for y in range(0, list_hbond_acceptors[z].shape[0]): # to go through the H-bond acceptors list
                if list_hbond_donors[i]['resSeq'][j] != list_hbond_acceptors[z]['resSeq'][y]:
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
os.system(f"plumed driver --plumed {inp_dir}/plumed_1.dat --ixtc {args_folded_trajectory} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_1_folded.out")
print("Running plumed_1.dat on the unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_1.dat --ixtc {args_unfolded_trajectory} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_1_unfolded.out")
os.chdir(inp_dir)

print("Reading COLVAR file for the first filter ...")
df = read_colvar('COLVAR_first_filter', '')
# drop H bonds with average length > 0.4nm
df = df.drop(df[(df['meanF'] > 0.4) & (df['meanU'] > 0.4) ].index).dropna()
# to remove the potential contacts within the same residues that don't display any difference between folded and unfolded
# drop H bonds with average difference between folded and unfolded < 0.01 nm
# TODO: This is delicate. Might drop things we want to keep. We should
# i) print what is discarded
# ii) probably turn it off for folded/unfolded sims. But this will leave in a lot of garbage, so further considerations are needed  
df = df.drop(df[(df['FminusU'].abs() < 0.01)].index).dropna()
df = df.sort_values(by=['labels'], ascending=True).reset_index()

###############################################################################################3
# Second filter to keep only the relevant contacts

# writes PLUMED script to obtain the angle information from the relevant H-bonds in df
print("Writing plumed_2.dat file to obtain angle information ...")
output = open('plumed_2.dat', 'w')
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
        
        # TODO: can be compressed with a function
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
        # TODO: can be compressed with a function
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
# TODO: heavy processing and generates two files. Might be intelligent to:
# i) Not have duplicates or, if this fails
# ii) load in memory the file and print it out on the same file to avoid duplicates and reduce garbage
lines_seen = set() # holds lines already seen
outfile = open('plumed_2_noduplicate.dat', "w")
for line in open('plumed_2.dat', "r"):
    if line not in lines_seen: # not a duplicate
        outfile.write(line)
        lines_seen.add(line)
outfile.close()

print("Running plumed_2.dat file on folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_2_noduplicate.dat --ixtc {args_folded_trajectory} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_2_folded.out")
print("Running plumed_2.dat file on unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_2_noduplicate.dat --ixtc {args_unfolded_trajectory} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_2_unfolded.out")
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

# TODO
# sometimes more than 1 H is picked up from aa like K, giving too much weight to that component of the CV. We need to clean them up
# they should naturally be consecutive
### HERE
#previous_soft_hbond = [-1,-1]
#for index, row in df_soft_sigF.iterrows():
#    # get atoms
#    print(index, row)
#    label_name = row['labels'].split("_")
#    atom_a = label_name[2].split("-")[0]
#    atom_b = label_name[3]
#    res_a = full_topology.atom(int(atom_a)).residue.resSeq
#    res_b = full_topology.atom(int(atom_b)).residue.resSeq
#    if [res_a, res_b] == previous_soft_hbond:
#        print("{res_a} and {res_b} have common soft hbonds")
#    else:
#        previous_soft_hbond = [res_a, res_b]
#    print(atom_a,atom_b)
#    print(res_a,res_b)

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

#print("Running plumed_3.dat on folded trajectory ...")
os.chdir(args_folded_dir)
print("Running plumed_3.dat on folded trajectory ...")
os.system(f"plumed driver --plumed {inp_dir}/plumed_3.dat --ixtc {args_folded_trajectory} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_3_folded.out")
print("Running plumed_3.dat on unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_3.dat --ixtc {args_unfolded_trajectory} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_3_unfolded.out")
os.chdir(inp_dir)

#############################################################################
# Adding the side chains

# writes a PLUMED file for the side-chain contacts
print("Writing plumed_SC.dat for the side-chain contacts ...")
output = open('plumed_SC.dat', 'w')
output.write('\n')

# to calculate the geometric center of the side chains (SC)
for res in number_residues:
    p = full_topology.select("protein and residue " + str(res) + " and not backbone and not (type H)") + 1
    p = p.tolist()
    if len(p) !=0:
        s = "SC" + str(res) + ": CENTER ATOMS=" + str(p) + "MASS"
        output.write(s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","="))
        output.write('\n')
    else:
        p = full_topology.select("protein and residue " + str(res) + " and backbone and name CA") + 1
        s = "SC" + str(res) + ": CENTER ATOMS=" + str(p) + "MASS"
        output.write(s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","="))
        output.write('\n')
output.write('\n')

# CV of the distances between the SCs (adjacent SC are not taken into account)
if len(number_residues) > 2:
    for i in range(0, len(number_residues)):
        res = number_residues[i]
        for j in range(i+1, len(number_residues)):
            exclude_res = number_residues[j]
            # don't consider neighbouring aa
            if res + 1 == exclude_res: continue
            label = "contside" + str(res) + "-" + str(exclude_res) + ": "
            output.write(str(label)+"COORDINATION GROUPA=SC" + str(res) + " GROUPB=SC" + str(exclude_res) + "  SWITCH={RATIONAL D_0=0.0 R_0=0.80 NN=4 MM=8}")
            output.write('\n')
    output.write('\n')

# calculate the distance between the two center of masse of residues pairs
if len(number_residues) > 2:
    for i in range(0, len(number_residues)):
        res = number_residues[i]
        for j in range(i+1, len(number_residues)):
            exclude_res = number_residues[j]
            # don't consider neighbouring aa
            if res + 1 == exclude_res: continue
            label = "SC" + str(res) + "-" + str(exclude_res) + ": "
            output.write(str(label) + "CENTER ATOMS=SC" + str(res) + ",SC" + str(exclude_res))
            output.write('\n')
    output.write('\n')

output.write('PRINT STRIDE='+str(args_stride)+ ' ARG=* FILE=COLVAR_SC')
output.write('\n')
output.write('\n')
output.close()

print("Running plumed_SC.dat on the folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_SC.dat --ixtc {args_folded_trajectory} --mc {args_mcfile} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_SC_folded.out")
print("Running plumed_SC.dat on the unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_SC.dat --ixtc {args_unfolded_trajectory} --mc {args_mcfile} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_SC_unfolded.out")
os.chdir(inp_dir)

# read the COLVAR files for the side chain contacts
print("Reading the COLVAR file for the side chains ...")
df_sc = read_colvar('COLVAR_SC', 'contside')
dffilteredFmU = df_sc[df_sc['FminusU'] > 0.0].dropna()
dffilteredFmU = dffilteredFmU[dffilteredFmU['meanF'] > args_cutoff].dropna()
dffilteredFmU = dffilteredFmU[dffilteredFmU['lda'] > args_lda].dropna()
dffilteredUmF = df_sc[df_sc['FminusU'] < 0.0].dropna()
dffilteredUmF = dffilteredUmF[dffilteredUmF['meanU'] > args_cutoff].dropna()
dffilteredUmF = dffilteredUmF[dffilteredUmF['lda'] > args_lda].dropna()
dffiltered = pd.concat([dffilteredFmU, dffilteredUmF])

##################################################################################
# Adding the solvation bias

# this PLUMED script will compute the solvation for each oxygen, nitrogen and carbon atoms

print("Writing plumed_solvation.dat for the solvation features ...")
output = open('plumed_solvation.dat', 'w')
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
os.system(f"plumed driver --plumed {inp_dir}/plumed_solvation.dat --ixtc {args_folded_trajectory} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_solvation_folded.out")
print("Running plumed_solvation.dat on the unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_solvation.dat --ixtc {args_unfolded_trajectory} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_solvation_unfolded.out")
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
print("Writing the final plumed file...")
output = open('plumed_final.dat', "w")
output.write('# This PLUMED file has been generated by bioinspired_features_generation.py on '+ str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'.\n')
#output.write('\n')
output.write('# Command : ' + full_command + '\n\n')
#output.write('\n')
#output.write('\n')
output.write("#RESTART\n")
#output.write('\n')
output.write(f'LOAD FILE={inp_dir}/LinearGhost.so\n')
#output.write('\n')
output.write(f'LOAD FILE={inp_dir}/CoordinationMapOMP.so\n')
#output.write('\n')
output.write(f'LOAD FILE={inp_dir}/HBondsOMP.so\n')
#output.write('\n')
#output.write('\n')
output.write('WO: GROUP ATOMS='+str(first_oxygen)+'-'+str(last_oxygen)+':'+str(water_model))
#output.write('\n')
output.write('WH: GROUP ATOMS='+str(first_oxygen+1)+'-'+str(last_oxygen+1)+':'+str(water_model)+','+str(first_oxygen+2)+'-'+str(last_oxygen+2)+':'+str(water_model)+'\n')
#output.write('\n')
output.write(f'rmsd_ca: RMSD REFERENCE={inp_dir}/{outname_protein_CA_ref} TYPE=OPTIMAL\n\n')
#output.write('\n')
#output.write('\n')

output.write("#############################################################\n")
#output.write('\n')
output.write("### Hydrogen Bonds ##########################################\n")
#output.write('\n')
output.write("#############################################################\n\n")
#output.write('\n')
#output.write('\n')

# TODO
# Probably change the local definition of purpose that is dangerous and pass it directly in the function call

# purpose is used in the process_HB function. If:
# 1 it prints LINEARGHOS
# 2 it prints HBONDS
# 3 it prints the water coordination for HBF VH-WO
# 4 it prints the water coordination for HBF VA-WH 
# 5 it prints the water coordination for HBU VH-WO
# 6 it prints the water coordination for HBU VA-WH

purpose = 1
process_HB(df_hard_soft,df_hard_sigF,df_hard_sigU,df_soft_sigF,df_soft_sigU,output,purpose)

#HBONDS
purpose = 2
#output.write('\n')
output.write('HB_sum: HBONDS...\n')
#output.write('\n')
process_HB(df_hard_soft,df_hard_sigF,df_hard_sigU,df_soft_sigF,df_soft_sigU,output,purpose)
output.write('   SUM\n')
#output.write('\n')
output.write('...\n')
#output.write('\n')

#COORDINATIONMAP
if df_hard_sigF.empty == False or df_soft_sigF.empty == False:
    #output.write('\n')
    output.write('\nNNC_HB_F: COORDINATIONMAP ...\n')
    #output.write('\n')
    process_virtual(df_virtualF, [df_hard_sigF,df_soft_sigF], output, prefix_HD="NPA", prefix_AB="NPD") 
    output.write('   SUM\n')
    #output.write('\n')
    output.write('...\n')
    #output.write('\n')
if df_hard_sigU.empty == False or df_soft_sigU.empty == False:
    #output.write('\n')
    output.write('\nNNC_HB_U: COORDINATIONMAP ...\n')
    #output.write('\n')
    process_virtual(df_virtualU, [df_hard_sigU,df_soft_sigU], output, prefix_HD="NPA", prefix_AB="NPD") 
    output.write('   SUM\n')
    #output.write('\n')
    output.write('...\n')
    #output.write('\n')

#Water coordination
purpose = 3
if df_hard_sigF.empty == False:
    #output.write('\n')
    output.write('\nNWO_F: COORDINATION GROUPA=')
    process_HB(df_hard_soft,df_hard_sigF,df_hard_sigU,df_soft_sigF,df_soft_sigU,output,purpose)
    purpose = 4
    #output.write('\n')
    output.write('\nNWH_F: COORDINATION GROUPA=')
    process_HB(df_hard_soft,df_hard_sigF,df_hard_sigU,df_soft_sigF,df_soft_sigU,output,purpose)
if df_hard_sigU.empty == False:
    purpose = 5
    #output.write('\n')
    output.write('\nNWO_U: COORDINATION GROUPA=')
    process_HB(df_hard_soft,df_hard_sigF,df_hard_sigU,df_soft_sigF,df_soft_sigU,output,purpose)
    purpose = 6
    #output.write('\n')
    output.write('\nNWH_U: COORDINATION GROUPA=')
    process_HB(df_hard_soft,df_hard_sigF,df_hard_sigU,df_soft_sigF,df_soft_sigU,output,purpose)

#HB CV, assuming that there are some HB contacts, hence cmapH_sum is not empty
coeff = [1.0, -1.0, -0.125, -0.0625, 1.0 if args_symm else 0.0, 0.125, 0.0625]

#output.write('\n')
#output.write('\n')
output.write('\n\ndiffHB: COMBINE ARG=HB_sum,')
if df_hard_sigF.empty == False or df_soft_sigF.empty == False:
    output.write('NNC_HB_F,')
    if df_hard_sigF.empty == False:
        output.write('NWO_F,NWH_F,')
if df_hard_sigU.empty == False or df_soft_sigU.empty == False:
    output.write('NNC_HB_U,')
    if df_hard_sigU.empty == False:
        output.write('NWO_U,NWH_U,')
output.write(' COEFFICIENTS='+str(coeff[0])+",")
if df_hard_sigF.empty == False or df_soft_sigF.empty == False:
    output.write(str(coeff[1])+",")
    if df_hard_sigF.empty == False:
        output.write(str(coeff[2])+","+str(coeff[3])+",")
if df_hard_sigU.empty == False or df_soft_sigU.empty == False:
    output.write(str(coeff[4])+",")
    if df_hard_sigU.empty == False:
        output.write(str(coeff[5])+","+str(coeff[6])+",")

output.write(' PERIODIC=NO\n\n')
#output.write('\n')
#output.write('\n')

##########
#####SC 
##########
output.write("#############################################################\n")
#output.write('\n')
output.write("### Side Chains #############################################\n")
#output.write('\n')
output.write("#############################################################\n\n")
#output.write('\n')
#output.write('\n')

# to calculate the center of mass of the side chains (SC), all are usually needed
for res in number_residues:
    p = full_topology.select("protein and residue " + str(res) + " and not backbone and not (type H)") + 1
    p = p.tolist()
    if len(p) !=0: #normal residues with side chains
        s = "SC" + str(res) + ": CENTER ATOMS=" + str(p) + "MASS"
        output.write(s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","=") + '\n')
        #output.write('\n')
    else:
        p = full_topology.select("protein and residue " + str(res) + " and backbone and name CA") + 1
        s = "SC" + str(res) + ": CENTER ATOMS=" + str(p) + "MASS"
        output.write(s.replace("]"," ").replace("["," ").replace(", ",",").replace("= ","=") + '\n')
        #output.write('\n')
output.write('\n')

# TODO: put purpose in the call of the function
#purpose is used in the process_SC function. If:
#1 it prints the SC CENTER virtual atoms
#2 it prints the SC CONTACTMAP
#3 it prints the SC NNC F
#4 it prints the SC NNC U

#SC CENTER
purpose = 1
process_SC(dffilteredFmU,dffilteredUmF,output,number_residues,purpose)

#SC CONTACTMAP
purpose = 2
#output.write('\n')
output.write('\nSC_sum: CONTACTMAP ...\n')
#output.write('\n')
process_SC(dffilteredFmU,dffilteredUmF,output,number_residues,purpose)
output.write('   SUM\n')
#output.write('\n')
output.write('...\n')
#output.write('\n')

#COORDINATIONMAP
if dffilteredFmU.empty == False:
    purpose = 3
    #output.write('\n')
    output.write('\nNNC_SC_F: COORDINATIONMAP ...\n')
    #output.write('\n')
    process_SC(dffilteredFmU,dffilteredUmF,output,number_residues,purpose)
    output.write('   SUM\n')
    #output.write('\n')
    output.write('...\n')
    #output.write('\n')
if dffilteredUmF.empty == False:
    purpose = 4
    #output.write('\n')
    output.write('\nNNC_SC_U: COORDINATIONMAP ...\n')
    #output.write('\n')
    process_SC(dffilteredFmU,dffilteredUmF,output,number_residues,purpose)
    output.write('   SUM\n')
    #output.write('\n')
    output.write('...\n')
    #output.write('\n')

#output.write('\n')
output.write('\ndiffSC: COMBINE ARG=SC_sum,')
if dffilteredFmU.empty == False:
    output.write('NNC_SC_F,')
if dffilteredUmF.empty == False:
    output.write('NNC_SC_U')
output.write(' COEFFICIENTS=1.0,')
if dffilteredFmU.empty == False:
    output.write('-1.0,')
if dffilteredUmF.empty == False:
    output.write('1.0')

output.write(' PERIODIC=NO\n\n')
#output.write('\n')
#output.write('\n')

output.write("#############################################################\n")
#output.write('\n')
output.write("### Extra Hydration CVs #####################################\n")
#output.write('\n')
output.write("#############################################################\n\n")
#output.write('\n')
#output.write('\n')

#WATER COORDINATION for Extra CVs
list_atoms = []
for i in range(0,len(df_solvation_C)):
    label = df_solvation_C['labels'][i]
    serial = re.findall(r'(\d+)', label)
    output.write(str(label)+': COORDINATION GROUPA='+str(serial[0])+ ' GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=10 D_MAX=1.0} NLIST NL_CUTOFF=1.5 NL_STRIDE=20 \n')
    #output.write('\n')
    list_atoms.append(str(label))
for i in range(0,len(df_solvation_NO)):
    label = df_solvation_NO['labels'][i]
    serial = re.findall(r'(\d+)', label)
    output.write(str(label)+': COORDINATION GROUPA='+str(serial[0])+ ' GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.3 NN=6 MM=10 D_MAX=1.0} NLIST NL_CUTOFF=1.5 NL_STRIDE=20 \n')
    #output.write('\n')
    list_atoms.append(str(label))
output.write('\n\n')
#output.write('\n')


output.write('PRINT STRIDE=10 ARG=diffHB,diffSC FILE=COLVAR_diff\n\n')
#output.write('\n')
#output.write('\n')

output.close()
##############################################################

# TODO: remove sed calls, bad practice
print("Running the final plumed file on the folded trajectory ...")
os.chdir(args_folded_dir)
os.system(f'sed -i "s/NLIST/#NLIST/g" {inp_dir}/plumed_final.dat')
os.system(f"plumed driver --plumed {inp_dir}/plumed_final.dat --ixtc {args_folded_trajectory} --mc {args_mcfile} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_final_folded.out")
print("Running the final plumed file on the unfolded trajectory ...")
os.chdir(args_unfolded_dir)
os.system(f"plumed driver --plumed {inp_dir}/plumed_final.dat --ixtc {args_unfolded_trajectory} --mc {args_mcfile} --pdb {inp_dir}/{outname_protein_ref} 1> plumed_final_unfolded.out")
os.chdir(inp_dir)
os.system(f'sed -i "s/#NLIST/NLIST/g" {inp_dir}/plumed_final.dat')

#### Wrap up
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

print("\n")
print("Analysis finished.")
print('Number of hard H-bonds significant in folded state: '+ str(len(df_hard_sigF)))
print('Number of hard H-bonds significant in unfolded state: '+ str(len(df_hard_sigU)))
print('Number of soft H-bonds significant in folded state: '+ str(len(df_soft_sigF)))
print('Number of soft H-bonds significant in unfolded state: '+ str(len(df_soft_sigU)))
print('Number of side chain contacts significant in folded state: '+ str(len(dffilteredFmU)))
print('Number of side chain contacts significant in unfolded state: '+ str(len(dffilteredUmF)))
print("\n")

# output pymol session if requested
if args_pymol:
    # load the protein and show it as licorice
    cmd.load(args_reference_protein)
    cmd.show_as("cartoon", "all")
    licorice_res = "residue " + str(number_residues[0])
    for i in range(1, len(number_residues)): licorice_res += "," + str(number_residues[i])
    cmd.show("licorice", licorice_res) # show only subsection in licorice
    cmd.hide("(all and hydro and (elem C extend 1))") # hide non-polar hydrogens
    cmd.hide("(solvent and (All))") # hide waters

    # loop on Hard-Folded, Hard-Unfolded, Soft-Folded, Soft-Unfolded hbonds
    for index, row in df_hard_soft.iterrows():
        # get atoms
        label_name = row['labels'].split("_")
        atom_a = label_name[2].split("-")[0]
        atom_b = label_name[3]
        # get type of hbond and prepare label name
        strength = "Hard" if row['Hard'] else "Soft"
        state = "Folded" if row ['Significant folded'] else "Unfolded"
        dist_name = f"hbond_{strength}_{state}_{atom_a}_{atom_b}" + "_lda_{0:.2f}".format(row['lda'])
        # draw the distance, hide the label (confusing), and color accordingly
        cmd.distance(dist_name, f"id {atom_a} and not solvent", f"id {atom_b} and not solvent")
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
