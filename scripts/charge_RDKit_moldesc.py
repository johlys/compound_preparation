#!/usr/bin/env python

#file with smiles, file with smirks, out_file
#range = 0 -> only one form per compound
#does not check stereo chemistry -> is not defined in Chemdiv files -> best to make the guesses after charging and tautomerizing the molecules
#if 3D input files -> stereo chemistry is kept
#if unique smile => has to be kekulised (no, that is not anylonger true, I fixed the "unique smiles" converter)
#this scripts also removes smaller fragments in complexes ("desalt")


#to do:
#1) test if chirality is properly assigned based on 3D coordinates
#2) test if SDfile is read properly
#4) Does the aromaticity model to be changed/adopted?
#5) are rings detected properly?
#6) check charges


from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdMolDescriptors
#import chemistry3 as chemistry # <--- this is only relevant for mol2-files, so we can skip it
#import my_mysql3 as mysql # <--- we will not use this
import sys, string
import pandas as pd
from argparse import ArgumentParser


#-----------

#*******************************************************************************
# MAIN help
#**********
#


description = "Script to charge molecules"
#usage = "drugpred.py [options]"

parser = ArgumentParser(description=description)

#-i works, --input_file does not work
parser.add_argument("-i", "--input_file", type=str, action="store", dest="infile", help="input file (SMILES (.smi) only)",required=True)
parser.add_argument("-o", "--output_file", type=str, action="store", dest="output_file",  help="output file (smiles)",required=True)
parser.add_argument("-c", "--charge", action="store_true", dest="charge",  default=False, help="write total charge to output file",required=False)
parser.add_argument("-p", "--ph", type=float, action="store", dest="ph", default= 7.0, help="pH, default 7.0",required=False)
parser.add_argument("-r", "--range", type=float, action="store", dest="ph_range", help="pH range",required=True)
parser.add_argument("-s", "--code",  action="store", dest="code", help="Name of field with supplier ID",required=False)
parser.add_argument("-d", "--debug",  action="store_true", dest="debug", default=False, help="Debug, default none",required=False)


args = parser.parse_args()
args_dict = vars(args)

print(args_dict)
locals().update(args_dict) #generates local variables from key / value pairs

infile_list = infile.split('.')
file_type = infile_list[-1]
print("type", file_type)
if file_type == 'smi' or file_type == 'csv':
	mol_list = Chem.SmilesMolSupplier(infile, delimiter="\t")
elif file_type == 'mol2': #this does not seem to be a well supported format
	print("reading mol2, currently not supported")
else:
	print ('unknown file type')
	sys.exit()

#rules_file = open(sys.argv[2], 'r')
out_file = open(output_file, 'w')


#print (charge, ph, ph_range)


#--------------------------------------------
#TO HANDLE TIMEOUT
def handler(signum, frame):
	print ("Timeout!")
	raise Exception("Timeout!")

#--------------------------------------------
def get_rules(): # <-- get this from definition_prot_rules.csv instead of MySQL
	prot_rules_github_url = "https://raw.githubusercontent.com/kjemist/compound_preparation/master/definition_prot_rules.csv"
	df_rules = pd.read_csv(prot_rules_github_url, index_col=0)
	rules = []
	for idx, row in df_rules.iterrows():
		rules.append([[row["smirks"]],row["type"],row["pka"],row["comment"], row.index])

	return rules
#--------------------------------------------
def reformat_rules(rules):
#OPENEYE has got a problem with writting [NH3+], ...
	rule_counter = 0
	for rule in rules:
		trans_counter = 0
		for transformation in rule[0]:
			educt, product = transformation.split('>>')
			#primary amines
			product = product.replace('[NH3+]', '[N+]([H])([H])[H]')
			#secundary amines
			product = product.replace('[NH2+]', '[N+]([H])([H])')
			#tertary amines
			product = product.replace('[NH+]', '[N+]([H])')
			product = product.replace('[nH+]', '[n+]([H])')
			product = product.replace('[nH]', '[n]([H])')
			reaction = educt + '>>' + product
			#print reaction
			rules[rule_counter][0][trans_counter] = reaction
			trans_counter = trans_counter + 1
		rule_counter = rule_counter + 1
	return rules
#--------------------------------------------

#when the molecule can react a lot, the script can hang => prepare time out
#signal.signal(signal.SIGALRM, handler)


#conn=mysql.connect2server('', 'webuser','purchasable')  	

#mol = OEGraphMol()


reaction = 0
type = 1
pka =2
comment = 3

#rules = read_rules(rules_file)
rules = get_rules()

#check if forward or backward reactions
counter = 0
for rule in rules:
	#print rule
	if (rule[type] == 'acid'):
		#if range =0 -> only one reaction
		if (rule[pka] < (ph - ph_range)) or (ph_range == 0 and rule[pka] <= (ph - ph_range)):
			direction = 'forward'
		elif rule[pka] > (ph + ph_range):
			direction = 'backwards'
		else:
			direction = 'both'
	else:
		if (rule[pka] > (ph + ph_range)) :
			direction = 'forward'
		elif rule[pka] < (ph - ph_range) or (ph_range == 0 and rule[pka] >= (ph + ph_range)):
			direction = 'backwards'
		else:
			direction = 'both'

	#print direction
	if direction != 'forward':
		#print rule[reaction][0]
		start = rule[reaction][0].find('>')
		backward_direction = rule[reaction][0][start+2:] +  '>>' + rule[reaction][0][:start]
	if direction == 'backwards':
			rules[counter][reaction][0] = backward_direction
	elif direction == 'both':
			rules[counter][reaction].append(backward_direction)
		#print 'test', rules[counter][reaction]
	counter = counter + 1
rules = reformat_rules(rules)
print("REFORMAT_RULES SUCCESFUL")
#print (rules)

print(mol_list)
for mol in mol_list:
	print(mol)
	print(mol.GetNumAtoms())
	print(rdMolDescriptors.CalcMolFormula(mol))
	print("TEST1")
#while OEReadMolecule(ifs, mol):
  #print mol.GetTitle()

	if file_type == 'sdf':
		print("Getting properties from sdf")
		try:
			name = mol.GetProp(code)
			name = name.replace(' ', '')
			mol.SetProp("_Name",name)
		except:
			print ('something wrong with this molecule, check also if code for SDF was correct')
			#sys.exit()
			continue
		#print (code, name)

	#print (mol.GetProp("_Name"), 'name')
	try:
		frag_mol = Chem.GetMolFrags(mol,True,True)
	except:
		print ("Could not read molecule")
		continue
	if len(frag_mol) > 1: #molecule was fragmented, keep largest fragment, could also be done with salt removal and SMARTS list
		#print (Chem.MolToSmiles(mol), ' < old')
		length = 0
		keep_frag = frag_mol[0]
		for frag in frag_mol:
			if frag.GetNumAtoms() > length:
				length = frag.GetNumAtoms()
				keep_frag = frag
		name = mol.GetProp("_Name")
		mol = keep_frag
		mol.SetProp("_Name",name)
	#print (Chem.MolToSmiles(mol))
	Chem.AssignStereochemistry(mol)
	#print (Chem.MolToSmiles(mol))


	#smiles_list = [chemistry.CanSmi(mol,iso,kek,chiral)]
	smiles_list = [Chem.MolToSmiles(mol)]
	#print smiles_list
	for rule in rules:
		#print rule
		tmp_list = []
		#tmp_list_kek = []	#only for kekulized simles unique thing works
		for smiles in smiles_list:
			#apply rules
			#print smiles
			old_smi = smiles
			for transformation in rule[0]:
			#make mol
				work_mol = Chem.MolFromSmiles(smiles) #this makes also explict Hs implicit
				try:
					work_mol = Chem.rdmolops.AddHs(work_mol) #make hydrogen atoms explicit, otherwise OE will complain for rules with explicit hydrogen atoms
				except:
					print ("could not add hydrogen atoms, skipping ", name)
					break
				#reacted = True
				to_do_list = [work_mol]
				while len(to_do_list) > 0: #make sure that same reaction is carried out on input smiles as often as possible, is an issue if two times the same group is contained in the lig => is not deprotanted at the same time
					new_to_do = []
					new_to_do_smiles = []
					for work_mol in to_do_list:
						umr = rdChemReactions.ReactionFromSmarts(transformation)
						try:
							product = umr.RunReactant(work_mol,0)
						except:
							print (smiles, mol.GetProp("_Name"), 'failed')
							print (transformation, rule[comment])
							work_mol = Chem.MolFromSmiles(smiles)
							print ('continue', smiles)
							continue #will need to set reacted = False to really go on

						counter = 0
						if len(product) == 0 and smiles not in tmp_list: #nothing has reacted, keep input smiles
							tmp_list.append(smiles)
						for ps in product: #which output does RunReactant acctually produce?
							for prod in ps:
								counter = counter + 1
								#print ('product: ', counter)
								#if 1 ==1:
								try:
									#print (smiles, mol.GetProp("_Name"))
									#print (Chem.MolToSmiles(prod))
									Chem.SanitizeMol(prod)	#check, if reaction has worked
								except: #transformatin has failed
								#if 2 ==2:
									print  ('----> error <------')
									print (smiles)
									print ('--> ' + mol.GetProp("_Name")  +', '+ transformation +', ' + rule[comment])
									#print (Chem.MolToSmiles(prod, sanitize=False))
									tmp_mol = prod
									tmp_mol = Chem.rdmolops.RemoveHs(tmp_mol,sanitize=False)
									tmp_smi = Chem.MolToSmiles(tmp_mol)
									print (tmp_smi, counter, product)
									print ('----> end <------')
									print
									if smiles not in tmp_list:
										tmp_list.append(smiles)	#save smiles
									continue
									#sys.exit()
								if prod.GetNumAtoms() + 3 < work_mol.GetNumAtoms(): #check if atoms got lost, one can happen, but be fixed
									print ("lost some atoms", mol.GetProp("_Name"))
									print ('--> ' + mol.GetProp("_Name")  +', '+ transformation +', ' + rule[comment])
									if smiles not in tmp_list:
										tmp_list.append(smiles)	#save smiles
									continue
								if umr.IsMoleculeReactant(prod): #check if molecule can react again => keep
									#print ("react again, *********")
									#convert to smiles to make comparision easyier
									prod = Chem.rdmolops.RemoveHs(prod)
									tmp_smi = (Chem.MolToSmiles(prod))								
									if tmp_smi not in new_to_do_smiles: #otherwise I have already got this, no need to save it twice
										new_to_do_smiles.append(tmp_smi)
										new_to_do.append(prod)
 
								#print (tmp_smi, reacted, counter)
									if debug:
										print ('----> reaction, can still react <------')	
										print (smiles)
										print ('--> ' + mol.GetProp("_Name")  +', '+ transformation +', ' + rule[comment])
										print ('Product:', tmp_smi)
										print ('----> end <------')
										print
								else: #can not reakt again
									prod = Chem.rdmolops.RemoveHs(prod) #convert to implict Hs
									tmp_smi = (Chem.MolToSmiles(prod)) #convert to implict Hs
							 
								#print (tmp_smi, reacted, counter)
									if tmp_smi not in tmp_list:
										tmp_list.append(tmp_smi)
										if debug:
											print ('----> reaction <------')	
											print (smiles)
											print ('--> ' + mol.GetProp("_Name")  +', '+ transformation +', ' + rule[comment])
											print ('Product:', tmp_smi)
											print ('----> end <------')
					to_do_list = new_to_do  #all single reactions for first molecules are carried out, now go on with the ones that can react again														   


		#update smiles_list
		smiles_list = []
		#print (len(tmp_list), 'len tmp_list')
		for i in tmp_list:
			smiles_list.append(i)
			#print (i), 'final'

	if debug:
		print ('results')
		print (smiles_list)

	print("SMILES_LIST")
	print(smiles_list)
	for moli in smiles_list:
		mol_out = Chem.MolFromSmiles(moli)
		Chem.SanitizeMol(mol_out) # this function does not return a new mol, but modifies the input mol
		smi_not_kek = Chem.MolToSmiles(mol_out)
		print(smi_not_kek)
		if charge == 'True':
				#create a molecule
				work_mol = OEGraphMol()
				OEParseSmiles(work_mol, moli)
				OEAssignImplicitHydrogens(work_mol)
				OEFormalPartialCharges(work_mol)
				total_charge = 0.0
				for atom in work_mol.GetAtoms():
						   total_charge = total_charge + atom.GetPartialCharge()
				if debug:
						print ('charge: ' + str(total_charge))
				out_file.write(smi_not_kek + '\t' + mol.GetProp("_Name") + '\t' + str(total_charge) + '\n')
		else:
				out_file.write(smi_not_kek + '\t' + mol.GetProp("_Name") + '\n')
	


		
print ('finished charger')

