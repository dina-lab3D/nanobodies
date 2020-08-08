from modeller import *
from modeller.automodel import *
from modeller import soap_loop
from modeller.scripts import complete_pdb

import Bio
from Bio import SeqIO
import cdr_annotation

import subprocess
import os,sys

current_home = os.path.dirname(sys.argv[0])
print (current_home)

blast_home = "/cs/labs/dina/dina/software/ncbi-blast-2.8.1+/bin/";

# runs blast
def run_blast(filename):
    out_file = filename + ".blast"
    blast_file = blast_home + "../pdbaa ";
    cmd = blast_home + "blastp -db " + blast_file + " -query " + filename + " >& " + out_file
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



# modeling script
if(len(sys.argv) != 2) :
    print ("Please provide input fasta_file")
    sys.exit(0)

fasta_file_name = sys.argv[1]
nb_sequence = get_sequence(fasta_file_name)
nb_sequence_id = get_sequence_id(fasta_file_name)
write_ali_file(nb_sequence)
#[cdr1_start, cdr1_end] = cdr_annotation.find_cdr1(nb_sequence)
[cdr2_start, cdr2_end] = cdr_annotation.find_cdr2(nb_sequence)
[cdr3_start, cdr3_end] = cdr_annotation.find_cdr3(nb_sequence)

# run blast and save top10 hits as templates
run_blast(fasta_file_name)
blast_file_name = fasta_file_name + ".blast"
count = 0
template_list = []

for line in open(blast_file_name):
    if "Chain" in line and count <10:
        fields = line.split(' ')
        fields2 = fields[0].split('_')
        if(len(fields2[1]) > 1):
            continue
        template = fields2[0].lower() + fields2[1]
        # ignore identical template
        if template in nb_sequence_id:
            continue

        template_list.append(template)

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

print (template_list)

# align templates using modeller salign function
log.verbose()
env = environ()
env.io.atom_files_directory = '.'
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

aln = alignment(env)
for (template_code) in (template_list):
    code = template_code[:4]
    chain = template_code[-1:]
    print (code,chain)
    mdl = model(env, file=template_code, model_segment=('FIRST:'+chain, '+120:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

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

# build 10 models with modeller automodel
model_num = 10
a = automodel(env, alnfile='nano-mult.ali',
              knowns=template_list,
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
model_num = 100 # can generate more if needed

# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        return selection(self.residue_range(cdr3_start+1, cdr3_end-2)) # focus

m = MyLoop(env,
           inimodel=best_model, # initial model of the target
           sequence='NANO')          # code of the target

m.loop.starting_model= 1           # index of the first loop model
m.loop.ending_model  = model_num   # index of the last loop model
m.loop.md_level = refine.slow      # loop refinement method; this yields
                                   # models quickly but of low quality;
                                   # use refine.slow for better models
m.make()

# score models
for n in range(1, 11):
    model_name = 'NANO.B999900%02d.pdb' % n
    mdl = complete_pdb(env, model_name)
    s = selection(mdl)   # all atom selection
    dope_score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',
                               normalize_profile=True, smoothing_window=15)
    soap_score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT',
                          normalize_profile=True, smoothing_window=15)
    rmsd = 0.0
    if os.path.isfile("ref.pdb"):
        cmd = "rmsd -t ref.pdb " + model_name + " | tail -n1 "
        rmsd_out = subprocess.check_output(cmd, shell=True)
        rmsd = float(rmsd_out.strip())
    print ("MODEL ", model_name, " dope-score: ", dope_score, " soap-score: ", soap_score, " rmsd: ", rmsd)
    f.write("MODEL "+ model_name + " dope-score: " + str(dope_score) + " soap-score: " + str(soap_score) + " rmsd: " + str(rmsd))

    # score loops
for i in range(1, model_num+1):
    # read model file
    code = "NANO.BL%04d0001.pdb" % i
    mdl = complete_pdb(env, code)
    s = selection(mdl)
    dope_score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT',
                               normalize_profile=True, smoothing_window=15)
    soap_score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT',
                          normalize_profile=True, smoothing_window=15)
    rmsd = 0.0
    if os.path.isfile("ref.pdb"):
        cmd = "rmsd -t ref.pdb " + code + " | tail -n1 "
        rmsd_out = subprocess.check_output(cmd, shell=True)
        rmsd = float(rmsd_out.strip())
    print ("LOOP ", code, " dope-score: ", dope_score, " soap-score: ", soap_score, " rmsd: ", rmsd)
    file_to_remove = "NANO.DL%04d0001" % i
    os.remove(file_to_remove)

# clean up templates
for (template_code) in (template_list):
    code = template_code[:4]
    chain = template_code[-1:]
    pdb = code + '.pdb'
    pdb_chain = template_code + '.pdb'
    os.remove(pdb)
    os.remove(pdb_chain)

# clean up other files
os.remove("NANO.rsr")
os.remove("NANO.ini")
os.remove("NANO.lrsr")
os.remove("NANO.sch")
os.remove("default")
