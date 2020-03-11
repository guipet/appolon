import numpy as np 
import nltk
import nltk.tree as Tree
from extraction import *
from PCYK import *
from tqdm import tqdm 

def flat_print(t):
	'''
	Input:
	- t: un arbre type nltk.tree
 	
 	Output:
 	- un str de l'arbre
	'''
	tag = t.label()
	if len(t) == 1:
		if type(t[0]) == str:
			if '+' in tag and tag != 'P+D':
				avant = tag[:tag.index("+")]
				après = tag[tag.index("+")+1 :]
				return '(' + avant + ' ' + '(' + après + ' ' + str(t[0]) + ')' + ')'
			else : 
				return '(' + tag + ' ' + str(t[0]) + ')'
		else:
			return '(' + tag + ' ' + flat_print(t[0]) + ')'
	else:
		s = []
		for i in range(len(t)):
			s.append(flat_print(t[i]))
		return '(' + tag + ' ' + ' '.join(s) + ')'


#on déparse le test
seperator = ' '
phrase_test = []
with open("test.txt","r") as file:

	for phrase in file:
			
		t = nltk.tree.Tree.fromstring(phrase)
		sent = seperator.join(t.leaves()) 
		phrase_test.append(sent)


#Initialisation de nos classes
grammer = PCFG()
grammer.build_pcfg('train.txt')
oov = OOV(grammer)

#on applique notre PCYK à ces phrases
test_deparse = open("evaluation_data.parser_output.txt","w")
for sent in tqdm(phrase_test):
	test = CYK(sent, grammer, oov)
	t = nltk.tree.Tree.fromstring(test)
	t.un_chomsky_normal_form(unaryChar='_')
	final = flat_print(t)
	print("la phrase finale est" + '\n')
	print(final)
	test_deparse.write( final + '\n')
test_deparse.close()









