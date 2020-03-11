import numpy as np 
from pcfg import *
from OOV import *
import nltk
import nltk.tree as Tree

def build_tree(words, back, i, j, A):
	'''
	Inputs:
	- words : list de str
	- back : dictionnaire backpointeurs
	- i : ième élément de la liste words
	- j : jième élément de la liste words
	-A : règle nltk
	'''

	if i == j:
		return( "".join([str(A), ' ', str(words[i])]))

	else : 
		B, C, s = back[i, j, A]
		return( "".join([str(A), ' (', build_tree(words, back, i, s, B), ') (', build_tree(words, back, s + 1, j, C), ')']) )




def CYK(sentence, grammer, oov):
	'''
	Inputs:
	- sentence : str 
	- grammer : classe de la grammaire
	- oov : class des oov

	Outputs:
	un str représentant la phrase parsée
	'''

	#on sépare en une liste de mot
	if type(sentence) != str:
		print("Erreur. Il faut que sentence soit un string")
	words = sentence.strip().split(' ')

	#variables
	table = {}
	back  = {}
	n = len(words)
	N = grammer.N

	#Initialisation
	for i in range(n):
		word = words[i]

		if not word in grammer.sigma:

			if i == 0 and i == (n-1) : #phrase à 1 mot
				prev = None
				nex  = None
			
			elif i == 0 and i != (n-1) : #début de phrase
				prev = None
				nex  = words[i+1]

			elif i == (n-1):  #fin de phrase
				prev = words[i-1]
				nex  = None

			else:
				prev = words[i-1]
				nex  = words[i+1] 

			print( "le mot original est", word)
			word = oov.closest_word(word, prev, nex)
			print("le remplaçant est ", word)

		for regle in N:
			if regle in grammer.lexicon.keys():

				if word in grammer.lexicon[regle].keys():
					table[i,i,regle] = grammer.lexicon[regle][word]

				else:
					table[i,i,regle] = 0

			else:
				table[i,i,regle] = 0

	#Récursion
	for l in range(1,n):
		for i in range(n-l) : 
			j = i+l

			for A in N:
				max_score = 0
				args = None

				if A in grammer.gram.keys():
					for couple in grammer.gram[A].keys():
						B, C = couple

						for s in range(i,j):

							if table[i, s, B] > 0 and table[s+1, j, C] > 0:
								score = table[i, s, B] * table[s+1, j, C] * grammer.gram[A][couple]

								if max_score < score : 
									max_score = score
									args = B, C, s

					if max_score > 0:
						back[i, j, A] = args

					table[i, j, A] = max_score

				else: 
					table[i,j,A] = 0
					back[i, j, A] = None


	##on essaye de refaire l'arbre
	max_score = 0
	args = None

	for A in N:
		if max_score < table[0, n-1, A]:
			max_score = table[0, n-1, A]
			args = 0, n-1, A

	if args == None:
		return( '(SENT (UNK))')

	else:

		return( '( SENT (' + build_tree(words, back, *args) + '))' )








 



