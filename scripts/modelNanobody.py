from modeller import *
from modeller.automodel import *
from modeller import soap_loop
from modeller.scripts import complete_pdb

import Bio
from Bio import PDB
from Bio import SeqIO
from Bio.Blast import NCBIXML
import cdr_annotation

import subprocess
import os,sys,getopt

current_home = os.path.dirname(sys.argv[0])
print (current_home)


blast_home = "/cs/labs/dina/dina/software/ncbi-blast-2.8.1+/bin/"
rmsd_prog = "/cs/staff/dina/utils/rmsd"
get_pdb = "/cs/staff/dina/scripts/getPDB.pl"
get_pdb_chains = "/cs/staff/dina/scripts/getPDBChains.pl"
renumber = "/cs/staff/dina/utils/srcs/renumber/renumber"
rmsd_align = "/cs/staff/dina/scripts/alignRMSD.pl"
get_frag_chain = "/cs/staff/dina/utils/get_frag_chain.Linux"

TEMPLATES_NUM = 10
FILTER_DIFF = 7.5
print(TEMPLATES_NUM)

# runs blast
def run_blast(filename):
    out_file = filename + ".blast"
    blast_file = blast_home + "../pdbaa "
    cmd = blast_home + "blastp -db " + blast_file + " -query " + filename + " -outfmt 5 >& " + out_file
    print (cmd)
    # Put stderr and stdout into pipes
    proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    return_code = proc.wait()
    if return_code != 0:
        print ("blast subprocess failed with exit code " , return_code)


# get first sequence
def get_sequence(fasta_filename):
    for seq_record in SeqIO.parse(fasta_filename, "fasta"):
        sequence = str(seq_record.seq)
        return sequence


# get first sequence id
def get_sequence_id(fasta_filename):
    for seq_record in SeqIO.parse(fasta_filename, "fasta"):
        sequence_id = str(seq_record.id)
        return sequence_id


# input for modeller instead of fasta
def write_ali_file(sequence):
    f = open("nano.ali", "w")
    f.write(">P1;NANO\n")
    f.write("sequence:NANO:::::::0.00: 0.00\n")
    f.write(sequence)
    f.write("*\n")
    f.close()


# get num_templates template pdbs based on blast xml
# returns a list of tuples (pdb_code, chain_id, template_start_res, template_end_res)
def get_templates_blast_xml(blast_file_name, num_templates = TEMPLATES_NUM, filter = False):
    result = open(blast_file_name,"r")
    records = NCBIXML.parse(result)
    item = next(records)
    counter = 0
    returned_templates = []
    seq_identity_list = []
    seq_file = open("seq_identity", 'w')

    for alignment in item.alignments:
        for hsp in alignment.hsps:
            if hsp.expect <0.01:
                pdb = alignment.title.split("|")[3].lower()
                chain = alignment.title.split("|")[4][0]
                seq_identity = 100*hsp.identities/hsp.align_length
                # get rid if numbering is off
                if hsp.sbjct_start > 20:
                    print("skipping " + pdb + chain + " " + str(seq_identity) + " " + str(hsp.sbjct_start))
                    continue

                # filter by 95% sequence identity
                if filter and seq_identity > 95.0:
                    print("skipping " + pdb + chain + " " + str(seq_identity))
                    continue

                if counter < num_templates:
                    returned_templates.append((pdb, chain, hsp.sbjct_start, hsp.sbjct_end))
                    seq_identity_list.append(seq_identity)

                    # get PDB files - all + chain only
                    pdb_name = pdb + ".pdb"
                    chain_pdb = pdb + chain + ".pdb"
                    if os.path.isfile(chain_pdb) :
                        print ("file exists", chain_pdb)
                    else: # getPDBChain
                        cmd = get_pdb_chains + " " + pdb + ":" + chain
                        print (cmd)
                        proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                        return_code = proc.wait()
                        if return_code != 0:
                            print ("getPDBChain.pl subprocess failed with exit code " , return_code)


                    if os.path.isfile(pdb_name) :
                        print ("file exists",pdb_name)
                    else: # getPDB
                        cmd = get_pdb + " " + pdb
                        print (cmd)
                        proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                        return_code = proc.wait()
                        if return_code != 0:
                            print ("getPDB.pl subprocess failed with exit code " , return_code)

                    seq_file.write(str(seq_identity) + "\n")
                    counter+=1
        if counter >= num_templates:
            break

    returned_templates = filter_templates(returned_templates, seq_identity_list)
    return returned_templates


def filter_templates(templates, seq_list):
    returned_templates = []
    max_identity = max(seq_list)
    seq_file = open("seq_identity_filter", 'w')
    for index in range(len(seq_list)):
        if (max_identity - seq_list[index]) <= FILTER_DIFF:
            returned_templates.append(templates[index])
            seq_file.write(str(seq_list[index]) + "\n")
    return returned_templates


# get template pdbs based on blast (old version)
def get_templates(blast_file_name, num_templates = TEMPLATES_NUM, filter = False):
    count = 0
    template_list = []

    for line in open(blast_file_name):
        if "Chain" in line and count <num_templates:
            fields = line.split(' ')
            fields2 = fields[0].split('_')
            if(len(fields2[1]) > 1):
                continue
            template = fields2[0].lower() + fields2[1]

            # ignore identical template
            if filter:
                score = int(line[70:77])
                if(score > 200):
                    print("skipping " + template + str(score))
                    continue

            template_list.append(template)

            # get PDB file
            pdb_name = template + ".pdb"
            if os.path.isfile(pdb_name) :
                print ("file exists",pdb_name)
            else: # getPDB
                template2 = fields2[0].lower() + ":" + fields2[1]
                cmd = "/cs/staff/dina/scripts/getPDBChains.pl " + template2
                print (cmd)
                proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                return_code = proc.wait()
                if return_code != 0:
                    print ("getPDBChains.pl subprocess failed with exit code " , return_code)

            # also get the whole PDB
            pdb_name = fields2[0].lower() + ".pdb"
            if os.path.isfile(pdb_name) :
                print ("file exists",pdb_name)
            else: # getPDB
                cmd = "/cs/staff/dina/scripts/getPDB.pl " + fields2[0].lower()
                print (cmd)
                proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
                return_code = proc.wait()
                if return_code != 0:
                    print ("getPDB.pl subprocess failed with exit code " , return_code)

            count+=1
    return template_list


def calc_dist(pdb_model, start, end, file):
    parser = Bio.PDB.PDBParser()
    model1 = parser.get_structure("1", "ref_renumber.pdb")[0]["H"]
    model2 = parser.get_structure("2", pdb_model)[0][' ']
    file.write("PDB: " + pdb_model)
    for shift in range(-3, 0):
        file.write(" " + str(shift) + ": " + str(model1[start+shift]["CA"] - model2[start+shift]["CA"]))
    file.write(" start: " + str(model1[start]["CA"] - model2[start]["CA"]))
    for shift in range(1, 4):
        file.write(" " + str(shift) + ": " + str(model1[end+shift]["CA"] - model2[end+shift]["CA"]))
    file.write(" end: " + str(model1[end]["CA"] - model2[end]["CA"]) + "\n")


# modeling script
help_message = "modelNanobody.py <inputFastaFile> [-t ignores sim. seq.] [-m model#] [-l loop num]"

# options defaults
filter_similar_sequences = False
model_num = 10
loop_model_num = 100

# parse options
#try:
opts, args = getopt.getopt(sys.argv[1:],'thm:l:', ['test', 'help','model_num=','loop_model_num='])
print (opts)
print (args)

for opt, arg in opts:
    if opt == '-h':
        print (help_message)
        sys.exit()
    if opt in ("-t", "--test"):
        filter_similar_sequences = True
    if opt in ("-m", "--model_num"):
        model_num = int(arg)
    if opt in ("-l", "--loop_model_num"):
        loop_model_num = int(arg)

# get fasta file name
if len (args) != 1:
    print (help_message)
    sys.exit()
fasta_file_name = args[0]
print ("Args: seq_file=" + fasta_file_name + " test_mode=" + str(filter_similar_sequences)
       + " model_num=" + str(model_num) + " loop_model_num=" + str(loop_model_num))

nb_sequence = get_sequence(fasta_file_name)
nb_sequence_id = get_sequence_id(fasta_file_name)
write_ali_file(nb_sequence)
#[cdr1_start, cdr1_end] = cdr_annotation.find_cdr1(nb_sequence) TODO- fix the cdr1 bug
[cdr2_start, cdr2_end] = cdr_annotation.find_cdr2(nb_sequence)
[cdr3_start, cdr3_end] = cdr_annotation.find_cdr3(nb_sequence)

# run blast and save top10 hits as templates
run_blast(fasta_file_name)
blast_file_name = fasta_file_name + ".blast"
template_list = get_templates_blast_xml(blast_file_name, filter = filter_similar_sequences)


# align templates using modeller salign function
log.verbose()
env = environ()
env.io.atom_files_directory = '.'
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

template_list2 = []
aln = alignment(env)
for template in template_list:
    print (template)
    pdb_code = template[0]
    chain = template[1]
    #template_start = template[2]
    #template_end = template[3]
    code = pdb_code+chain
    #start = str(template_start)+':' +chain
    #end = str(template_end)+':'+ chain
    try:
        mdl = model(env, file=code, model_segment=('FIRST:'+chain, '+120:'+chain))
        aln.append_model(mdl, atom_files=pdb_code, align_codes=pdb_code+chain)
        template_list2.append(code)
    except OSError:  # pdb file not found...
        continue


for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, #write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

aln.write(file='templates.ali', alignment_format='PIR')
aln_block = len(aln)

# Read aligned sequence(s):
aln.append(file='nano.ali', align_codes='NANO')

# Structure sensitive variable gap penalty sequence-sequence alignment:
aln.salign(output='', max_gap_length=10,
           gap_function=True,   # to use structure-dependent gap penalty
           alignment_type='PAIRWISE', align_block=aln_block,
           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=99,
           gap_penalties_1d=(-450, 0),
           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
           similarity_flag=True, local_alignment=True)

aln.write(file='nano-mult.ali', alignment_format='PIR')

# scoring
sp = soap_loop.Scorer()

print ("TL",template_list2)

# build 10 models with modeller automodel
a = automodel(env, alnfile='nano-mult.ali',
              knowns=template_list2,
              sequence='NANO',assess_methods=(assess.DOPE,sp))
a.starting_model = 1
a.ending_model = model_num
a.make()

# score the models
best_score = 0
best_model = "x"
for n in range(1, model_num+1):
    model_name = 'NANO.B999900%02d.pdb' % n
    mdl = complete_pdb(env, model_name)
    s = selection(mdl)   # all atom selection
    score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT', #file='TCR.B9999000' + str(n) + '.profile',
                     normalize_profile=True, smoothing_window=15)
    file_to_remove = 'NANO.V999900%02d' % n
    os.remove(file_to_remove)
    file_to_remove = 'NANO.D000000%02d' % n
    os.remove(file_to_remove)
    if(score < best_score):
        best_score = score
        best_model = model_name

print (best_model, best_score)

# loop remodeling
# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        return selection(self.residue_range(cdr3_start-1, cdr3_end+2)) # focus

m = MyLoop(env,
           inimodel=best_model, # initial model of the target
           sequence='NANO')          # code of the target

m.loop.starting_model= 1           # index of the first loop model
m.loop.ending_model  = loop_model_num   # index of the last loop model
m.loop.md_level = refine.slow      # loop refinement method; this yields
                                   # models quickly but of low quality;
                                   # use refine.slow for better models
m.make()


# renumber for calculating cdrs rmsd (same length of models...)
if os.path.isfile("ref.pdb"):
    subprocess.run(renumber + " ref.pdb > ref_renumber.pdb", shell=True)

    #subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(cdr1_start) + " " + str(cdr1_end) + " > ref_cdr1.pdb", shell=True)
    subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(cdr2_start) + " " + str(cdr2_end) + " > ref_cdr2.pdb", shell=True)
    subprocess.run(get_frag_chain + " ref_renumber.pdb H " + str(cdr3_start) + " " + str(cdr3_end) + " > ref_cdr3.pdb", shell=True)


# score models
f = open("scores.txt", "w")
cdr_dist = open("cdrs_dist", "w")

for n in range(1, model_num+1):
    model_name = 'NANO.B999900%02d.pdb' % n
    mdl = complete_pdb(env, model_name)
    s = selection(mdl)   # all atom selection

    dope_score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',
                               normalize_profile=True, smoothing_window=15)
    soap_score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT',
                          normalize_profile=True, smoothing_window=15)
    rmsd = 0.0
    cdr1_rmsd = 0.0
    cdr2_rmsd = 0.0
    cdr3_rmsd = 0.0
    if os.path.isfile("ref.pdb"):

        cmd = rmsd_prog + " -t ref.pdb " + model_name + " | tail -n1 "
        rmsd_out = subprocess.check_output(cmd, shell=True)
        rmsd = float(rmsd_out.strip())

        # cdr 1,2,3 rmsd

        subprocess.run(rmsd_align + " ref.pdb" + " " + model_name, shell=True)  # align to get rmsd of cdr without cheating...

        #  subprocess.run(get_frag_chain + " " + model_name.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr1_start) + " " + str(cdr1_end) + " > temp_cdr1.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + model_name.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr2_start) + " " + str(cdr2_end) + " > temp_cdr2.pdb", shell=True)
        subprocess.run(get_frag_chain + " " + model_name.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr3_start) + " " + str(cdr3_end) + " > temp_cdr3.pdb", shell=True)

        #  cdr1_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr1.pdb temp_cdr1.pdb | tail -n1 ", shell=True).strip())
        cdr2_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr2.pdb temp_cdr2.pdb | tail -n1 ", shell=True).strip())
        cdr3_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr3.pdb temp_cdr3.pdb | tail -n1 ", shell=True).strip())

        calc_dist(model_name.replace(".pdb", "_tr.pdb"), cdr3_start, cdr3_end, cdr_dist)
        os.remove(model_name.replace(".pdb", "_tr.pdb"))


    print ("MODEL ", model_name, " dope-score: ", dope_score, " soap-score: ", soap_score, " rmsd: ", rmsd, " cdr1-rmsd: ", 0, " cdr2-rmsd: ", cdr2_rmsd, " cdr3-rmsd: ", cdr3_rmsd)
    ""
    f.write("MODEL "+ model_name + " dope-score: " + str(dope_score) + " soap-score: " + str(soap_score) + " rmsd: " + str(rmsd) +
            " cdr1-rmsd: " + str(cdr1_rmsd) + " cdr2-rmsd: " + str(cdr2_rmsd) + " cdr3-rmsd: " + str(cdr3_rmsd) + "\n")

# score loops
for i in range(1, loop_model_num+1):
    # read model file
    code = "NANO.BL%04d0001.pdb" % i
    mdl = complete_pdb(env, code)
    s = selection(mdl)
    dope_score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',
                               normalize_profile=True, smoothing_window=15)
    soap_score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT',
                          normalize_profile=True, smoothing_window=15)
    rmsd = 0.0
    cdr1_rmsd = 0.0
    cdr2_rmsd = 0.0
    cdr3_rmsd = 0.0

    if os.path.isfile("ref.pdb"):
        cmd = rmsd_prog + " -t ref.pdb " + code + " | tail -n1 "
        rmsd_out = subprocess.check_output(cmd, shell=True)
        rmsd = float(rmsd_out.strip())

    # cdr 1,2,3 rmsd

    subprocess.run(rmsd_align + " ref.pdb" + " " + code, shell=True)  # align to get rmsd of cdr without cheating...

    #subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr1_start) + " " + str(cdr1_end) + " > temp_cdr1.pdb", shell=True)
    subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr2_start) + " " + str(cdr2_end) + " > temp_cdr2.pdb", shell=True)
    subprocess.run(get_frag_chain + " " + code.replace(".pdb", "_tr.pdb") + " ' ' " + str(cdr3_start) + " " + str(cdr3_end) + " > temp_cdr3.pdb", shell=True)

    #cdr1_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr1.pdb temp_cdr1.pdb | tail -n1 ", shell=True).strip())
    cdr2_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr2.pdb temp_cdr2.pdb | tail -n1 ", shell=True).strip())
    cdr3_rmsd = float(subprocess.check_output(rmsd_prog + " ref_cdr3.pdb temp_cdr3.pdb | tail -n1 ", shell=True).strip())

    print ("LOOP ", code, " dope-score: ", dope_score, " soap-score: ", soap_score, " rmsd: ", rmsd, " cdr1-rmsd: ", 0, " cdr2-rmsd: ", cdr2_rmsd, " cdr3-rmsd: ", cdr3_rmsd)
    f.write("LOOP "+ code + " dope-score: " + str(dope_score) + " soap-score: " + str(soap_score) + " rmsd: " + str(rmsd) +
            " cdr1-rmsd: " + str(cdr1_rmsd) + " cdr2-rmsd: " + str(cdr2_rmsd) + " cdr3-rmsd: " + str(cdr3_rmsd) + "\n")

    file_to_remove = "NANO.DL%04d0001" % i
    os.remove(file_to_remove)
    os.remove(code.replace(".pdb", "_tr.pdb"))

# clean cdrs files
if os.path.isfile("ref.pdb"):
    # os.remove("temp_cdr1.pdb")
    os.remove("temp_cdr2.pdb")
    os.remove("temp_cdr3.pdb")
    os.remove("ref_renumber.pdb")

# clean up templates
for template in template_list:
    code = template[0]
    chain = template[1]
    pdb = code + '.pdb'
    pdb_chain = code + chain + '.pdb'
    if os.path.exists(pdb):
        os.remove(pdb)
        os.remove(pdb_chain)

# clean up other files
os.remove("NANO.rsr")
os.remove("NANO.ini")
os.remove("NANO.lrsr")
os.remove("NANO.sch")
os.remove("default")
