from modeller import *
from modeller.automodel import *
from modeller import soap_protein_od
from modeller.scripts import complete_pdb

import subprocess
import os,sys

current_home = os.path.dirname(sys.argv[0])
print current_home;

blast_home = "/cs/labs/dina/dina/software/ncbi-blast-2.8.1+/bin/";
nano_pdbs_folder = '/vol/sci/bio/bio3d/dina/BioSystems/nanobodies/pdbs/'

def run_blast(filename):
    out_file = filename + ".blast"
    blast_file = blast_home + "../pdbaa ";
    cmd = blast_home + "blastp -db " + blast_file + " -query " + filename + " >& " + out_file
    print cmd
    # Put stderr and stdout into pipes
    proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    return_code = proc.wait()
    if return_code != 0:
        print "blast subprocess failed with exit code " , return_code

if(len(sys.argv) != 4) :
    print "Please provide input fasta_file cdr3_start cdr3_end";
    sys.exit(0)

fasta_file_name = sys.argv[1]
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
        #print fields2[0],":",fields2[1]
        template = fields2[0].lower() + fields2[1]
        template_list.append(template)

        pdb_name = nano_pdbs_folder + template + ".pdb"
        if os.path.isfile(pdb_name) :
            print "file exists",pdb_name
        else: # getPDB
            template2 = fields2[0].lower() + ":" + fields2[1]
            cmd = "/cs/staff/dina/scripts/getPDBChains.pl " + template2
            print cmd
            proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
            return_code = proc.wait()
            if return_code != 0:
                print "getPDBChains.pl subprocess failed with exit code " , return_code
            pdb = template + ".pdb"
            print pdb, pdb_name
            os.rename(pdb, pdb_name)

        # also get the whole PDB
        pdb_name = nano_pdbs_folder + fields2[0].lower() + ".pdb"
        if os.path.isfile(pdb_name) :
            print "file exists",pdb_name
        else: # getPDB
            cmd = "/cs/staff/dina/scripts/getPDB.pl " + fields2[0].lower()
            print cmd
            proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
            return_code = proc.wait()
            if return_code != 0:
                print "getPDB.pl subprocess failed with exit code " , return_code
            pdb = fields2[0].lower() + ".pdb"
            print pdb, pdb_name
            os.rename(pdb, pdb_name)

        count+=1


print template_list;

log.verbose()
env = environ()
env.io.atom_files_directory = nano_pdbs_folder
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

aln = alignment(env)
for (template_code) in (template_list):
    code = template_code[:4]
    chain = template_code[-1:]
    print code;
    print chain;
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

f = open('gnuplotfile', 'w')
f.write('plot ')

sp = soap_protein_od.Scorer()
#for (template_code) in (template_list):
#    code = template_code[:4]
#    chain = template_code[-1:]
#    mdl = complete_pdb(env, code + chain + '_fit.pdb')
#    s = selection(mdl)   # all atom selection
#    score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT', #file=code + '.profile',
#                     normalize_profile=True, smoothing_window=15)
    #f.write('\''+ code + '.profile\' u 1:42 w l lc rgb \'#61a0e2\',')


a = automodel(env, alnfile='nano-mult.ali',
              knowns=template_list,
              sequence='NANO',assess_methods=(assess.DOPE,sp))
a.starting_model = 1
a.ending_model = 10
a.make()

best_score = 0
best_model = "x"
for n in range(1, 9):
    model_name = 'NANO.B9999000' + str(n) + '.pdb'
    mdl = complete_pdb(env, model_name)
    s = selection(mdl)   # all atom selection
    score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT', #file='TCR.B9999000' + str(n) + '.profile',
                     normalize_profile=True, smoothing_window=15)
    if(score < best_score):
        best_score = score
        best_model = model_name
#    f.write('\'NANO.B9999000' + str(n) + '.profile\' u 1:42 w l lc rgb \'#e26261\',')


model_name = 'NANO.B99990010.pdb'
mdl = complete_pdb(env, model_name)
s = selection(mdl)   # all atom selection
score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT', #file='TCR.B99990010.profile',
                     normalize_profile=True, smoothing_window=15)
if(score < best_score):
    best_score = score
    best_model = model_name
#f.write('\'NANO.B99990010.profile\' u 1:42 w l lc rgb \'#e26261\'')


print best_model, best_score


# loop remodeling
env.io.atom_files_directory = '.'
cdr3_start = sys.argv[2]
cdr3_end = sys.argv[3]


# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        return selection(self.residue_range(cdr3_start, cdr3_end))

m = MyLoop(env,
           inimodel=best_model, # initial model of the target
           sequence='NANO')          # code of the target

m.loop.starting_model= 1           # index of the first loop model
m.loop.ending_model  = 99          # index of the last loop model
m.loop.md_level = refine.slow      # loop refinement method; this yields
                                   # models quickly but of low quality;
                                   # use refine.slow for better models

m.make()


for i in range(1, 99):
    # read model file
    code = "NANO.BL%04d0001.pdb" % i
    mdl = complete_pdb(env, code)
    s = selection(mdl)
    score = s.assess_dope(output='ENERGY_PROFILE NO_REPORT', #file='TvLDH.profile',
                          normalize_profile=True, smoothing_window=15)
    print "LOOP ", code, " dope-score: ", score
