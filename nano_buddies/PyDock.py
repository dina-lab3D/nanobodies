
import argparse
import os
import subprocess
pd_home = "/cs/staff/dina/projects2/PatchDock/"
imp_home = "/cs/labs/dina/dina/libs/imp_cluster/"

# max memory for each batch
MEMORY = "10000m"

# max time for each batch
TIME = "72:0:0"

SCRIPT_PATH = "/cs/labs/dina/tomer.cohen13/nanobodies/nano_buddies/nanobodies_script.py"

# the begining of the script for cluster
INTRO = "#!/bin/tcsh\n" \
        "#SBATCH --mem="+ MEMORY +"\n" \
        "#SBATCH -c1\n" \
        "#SBATCH --time=" + TIME + "\n"


def main(directory):

    os.chdir(directory)
    antigen = os.path.basename(directory)

    antigenPDB = antigen + ".pdb"

    # my $curr = cwd;
    # `cp $antigen .`;
    antigen_chains = subprocess.run("~dina/scripts/chainSelector.pl " + antigenPDB, shell=True,
                                    capture_output=True, universal_newlines=True).stdout.replace("H", "").replace(" ", "").replace("\n", "")

    subprocess.run("~dina/utils/getChain.Linux " + antigen_chains + " " + antigenPDB + " > " + "antigen.pdb", shell=True)

# with open("dock_script.sh", 'w') as script_file:
#
#         script_file.write(INTRO)
#         script_file.write("cd " + os.getcwd())
#
#     for (my $j=0; $j < 5; $j++) {
#         my $loopfile = "nb_loop_" . $j . ".pdb";
#
#     # set nanobody chain id to H
#     my $cmd = "/cs/staff/dina/scripts/chainger.pl $loopfile ' ' 'H'";
#     print OUT "$cmd\n";
#
#     # parameter file
#     $cmd = "$pd_home/buildParams.pl $loopfile $antigenPDB 4.0 AA";
#     print OUT "$cmd\n";
#     $cmd = "mv params.txt params$j.txt";
#     print OUT "$cmd\n";
#     $cmd = "cp mycdrs3 cdrs3";
#     print OUT "$cmd\n";
#     $cmd = "cp myframe frame";
#     print OUT "$cmd\n";
#
#     # patchdock
#     if(not (-e "docking$j.res" and (-s "docking$j.res" > 1000))) {
#     $cmd = "$pd_home/patch_dock.Linux params$j.txt docking$j.res";
#     print OUT "$cmd\n";
#     $cmd = "$pd_home/PatchDockOut2Trans.pl docking$j.res > trans$j";
#     print OUT "$cmd\n";
#     }
#
#     # soap score
#     if(not (-e "soap_score$j.res" and (-s "soap_score$j.res" > 200))) {
#     $cmd = "$imp_home/setup_environment.sh $imp_home/bin/soap_score $antigenPDB $loopfile trans$j -o soap_score$j.res";
#     print OUT "$cmd\n";
#     }
#     }
#
#     close OUT;
#
#     `sbatch --time=8:0:0 dscript.sh`;

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="directory path containing the pdb directories")
    args = parser.parse_args()
    main(args.directory)

