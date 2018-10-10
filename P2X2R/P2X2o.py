from modeller import *
from modeller.automodel import *
log.verbose()

class RSmodel (automodel):
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['A','B','C'], renumber_residues=[25,25,25])   
        # start residues numbering from 22 instead from 1
    def special_restraints(self,aln):
        s1 = selection(self.chains['A']).only_atom_types('CA')
        s2 = selection(self.chains['B']).only_atom_types('CA')
        s3 = selection(self.chains['C']).only_atom_types('CA')        
        self.restraints.symmetry.append(symmetry(s1, s2, 1.0))
        self.restraints.symmetry.append(symmetry(s1, s3, 1.0))
        self.restraints.symmetry.append(symmetry(s2, s3, 1.0))
    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)
 
env = environ() 			# create new MODELLER environment to build this model in
env.io.hetatm = True		# read heteroatom entries 
# env.io.water = True

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# Use modified automodel class
a = RSmodel(env,
			alnfile='P2X2_P2X3o.pir',  	# alignment filename
			knowns=('5svk_t'), 		    # codes of the template(s)
			sequence='P2X2o',		    # code of the target
            assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1	#index of the first model
a.ending_model = 5		#index of the last model (use 5 for testing 50 or 100 for production) 
a.make()				# do the actual modeling


# Get a list of all successfully built models from a.outputs
ok_models = filter(lambda x: x['failure'] is None, a.outputs)

# Rank the models by DOPE score
def numSort(a,b,key = 'DOPE score'):
    return cmp(a[key],b[key])
ok_models.sort(numSort)

# Print list of models ranked according to DOPE score
output = open('ranked_models.txt','w')
output.write('Model file name\t molpdf score\t DOPE score\n')
for m in ok_models:
    results = '%s\t' % m['name']
    results += '%.3f\t' % m['molpdf']
    results += '%.3f\n' %m['DOPE score']
    output.write(results)
output.close()
    





