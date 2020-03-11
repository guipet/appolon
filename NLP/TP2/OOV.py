import numpy as np 
import pickle
import nltk
import nltk.tree as Tree
from pcfg import *
import re

# Noramlize digits by replacing them with #
DIGITS = re.compile("[0-9]", re.UNICODE)

#cosine similarity
def cos_similarity(word1, word2):
	'''
	Inputs:
	-word1 : vecteur 
	-word2 : vecteur

	Output: 
	-score : réel

	calcule la cos similarity entre deux mots
	'''
	if (np.linalg.norm(word1) == 0) or (np.linalg.norm(word2) == 0):
		return(0)

	else:
		score = (word1.dot(word2)) /  (np.linalg.norm(word1)*np.linalg.norm(word2))
		return(score)



def Damereau_Levenshtein(word1, word2):
	'''
	Inputs :
	-word1 : str
	-word2 : str

	Output:
	-Levenshtein distance

	Calcule la distance de Levenshtein entre 2 mots.
	'''
	n, m = len(word1), len(word2)
	D = np.zeros( (n+1,m+1))
	#Initialisation
	for i in range(n+1):
		D[i,0] = i
	for j in range(m+1):
		D[0,j] = j

	#récursion
	for i in range(1, n+1):
		for j in range(1,m+1):
			cost = 0 + 2*(word1[i-1] != word2[j-1])

			D[i,j] = min( [D[i-1, j] + 1, D[i, j-1] + 1, D[i-1, j-1] + cost])

			if i > 1 and j > 1 and word1[i-1] == word2[j-2] and word1[i-2] == word2[j-1]:
				D[i,j] = min( [D[i,j], D[i-2, j-2] + cost])

	return(D[n,m])


class OOV(object):

	def __init__(self, grammer):
		'''
		Input:
		-grammar class

		Outputs:
		- words_embed : dic des mots du pickle
		- embeddings : matrice des vecteurs d'embedding du pickle
		- vocab : dic des indices des mots du pickle
		- voc_lexicon : liste de str du vocabulaire du train
		- mot2id : dic des indices des mots du lexicon de grammar
		- id2mot : dic des mots liés à leurs indices du lexicon de grammar
		- bigram_matrix : matrice du bigramme du lexicon
		- unigram : vecteur des proba de l'unigramme


		Renvoie le mot le plus proche du train pour remplacer le OOV
		'''

		self.grammer = grammer
		self.words_embed, self.embeddings = pickle.load(open('polyglot-fr.pkl', 'rb'), encoding='latin1')
		#vocabulaire des mots du polyglot
		self.vocab = {}
		for idx, w in enumerate(self.words_embed):
			self.vocab[w] = idx


		#dictionnaires correspondant au mot du lexicon
		self.mot2id = {}   #retrouver l'indice du mot
		self.id2mot = {}   #retrouver un mot selon son indice

		for i, mot in enumerate(self.grammer.sigma):
			self.mot2id[mot] = i
			self.id2mot[i] = mot

		self.mot2id['<s>'] = len(self.mot2id)
		self.id2mot[ len(self.grammer.sigma) ] = '<s>'


        #on crée la matrice de bigramme et l'unigramme
		print("création matrice bigramme...")
		n = len( self.grammer.sigma )
		total = 0
		self.bigram_matrix = np.ones ((n+1,n+1)) #Laplace smoothing

		with open('train.txt', "r") as txt:
			for phrase in txt:

				t = nltk.tree.Tree.fromstring(phrase)

				#on rajoute des balises de début et de fin (SI JAMAIS FAUTE EN DÉBUT OU FIN, CELA PEUT AMENER À DES RÉSULTATS ERRONÉS)
				feuilles = ['<s>']
				for mot in t.leaves():
					feuilles.append(mot)
				feuilles.append('<s>')

				#on uupdate le nombre de mot
				total += len(feuilles[1: len(feuilles)-1])
				#on update la matrice du bigramme 
				for idx in range(len(feuilles) - 1):

					ligne = feuilles[idx]
					col   = feuilles[idx + 1]
					self.bigram_matrix[ self.mot2id[ligne], self.mot2id[col] ] +=1

			normalisation = np.sum(self.bigram_matrix, axis = 1)
			self.bigram_matrix = (self.bigram_matrix.T / normalisation).T
			self.bigram_matrix = np.nan_to_num(self.bigram_matrix)
		print('matrice bigramme créée')

		#on crée l'unigramme
		print('création unigramme...')

		self.unigram = np.zeros(n)
		for mot in self.grammer.count_words.keys():
			self.unigram[ self.mot2id[mot] ] = grammer.count_words[mot] / total 

		print('unigramme crée')

    	#il faudrait créer l'unigramme et appliquer la formule des slides.
    	#comment bien paramétrer lambda


	def case_normalizer(self, word, dico):
		'''
    	Inputs:
    	-word : str qu'on veut modifier
    	- dico : soit vocab si on veut comparer aux mot de l'embedding soit mot2id si on veut comparer au lexique

    	Ouput:
    	- str proche de word s'il existe, sinon on retourne le mot 
    	'''
		w = word
		lower = (dico.get(w.lower(), 1e12), w.lower())
		upper = (dico.get(w.upper(), 1e12), w.upper())
		title = (dico.get(w.title(), 1e12), w.title())
		results = [lower, upper, title]
		results.sort()
		index, w = results[0]
		if index != 1e12:
			return(w)
		return(word)


	def normalize(self, word, dico):
		'''
    	Input:
    	-word : str qu'on veut normaliser
    	-dico : dict selon lequel on normalise
    	'''

    	#on remplace les chiffres par des #

		if not word in dico:
			word = DIGITS.sub("#", word)
    	#si toujours pas dans word, alors on modofie les polices
		if not word in dico:
			word = self.case_normalizer(word, dico)
		return(word)


	def candidats_max(self, word, k = 8):
		'''
		Inputs:
		- word : str dont on cherche les candidats
		- k : un réel déterminant la distance maximale pour Levenshtein
		'''
		candidats = []
		dist = []
		for mot in self.mot2id.keys():
			if mot != '<s>':  #peut arriver comme la distance est grande
				distance = Damereau_Levenshtein(word, mot)

				if distance <= k:
					candidats.append(mot)
					dist.append(distance)

		return(candidats, dist)

	def proba_bigram(self,prev, oov, nex):
		'''
    	Inputs:
    	- prev : str. Mot précedant le oov dans la phrase 
    	- oov : str. Mot inconnu
    	- next : str. Mot succédant le oov dans la phrase

    	Output:
    	- les probas correspondantes.
    	'''
		début = self.bigram_matrix.shape[0]-1
		fin   = self.bigram_matrix.shape[1]-1

		if prev == None:  #début de phrase
			prev_oov = self.bigram_matrix[début, self.mot2id[oov]]

		elif prev in self.mot2id.keys(): #si dans le lexicon, on prend compte le contexte
			prev_oov = self.bigram_matrix[self.mot2id[prev], self.mot2id[oov]]

		else: #sinon on ne sait pas l'impact du contexte
			prev_oov = 1

		if nex == None:
			oov_nex = self.bigram_matrix[self.mot2id[oov], fin]

		elif nex in self.mot2id.keys():
			oov_nex = self.bigram_matrix[self.mot2id[oov], self.mot2id[nex]]
		else:
			oov_nex = 1

		return( prev_oov * oov_nex) 


	def closest_word(self, word, prev, nex, k = 2, coef = 10000000):  #nos probas valent en générale 10**(-8)
		'''
    	Inputs:
    	- word : str qui n'est pas dans le lexique
    	- prev : str précédant word dans la phrase
    	- nex : str succédant word dans la phrase
    	- k : réel 
    	- coef : coef d'importance au modèle langage.

    	Output:
    	- str le plus proche dans le lexique
    	'''

    	#1ére étape, on normalise le mot
		w = self.normalize(word, self.mot2id) 

		#si la normalisation est dans le lexique
		if w in self.mot2id.keys():
			return(w)


    	#Si la normalisation est dans le polyglot + candidat du train si jamais erreur d'écriture
		w = self.normalize(word, self.vocab)
		if w in self.vocab.keys():

			#cosine similarity avec le train
			print("la normalisation est dans le polyglot. On verifie aussi les erreurs d'écriture")
			embedding = self.embeddings[ self.vocab[w] ]  #on prend l'embedding du mot dans le polyglot
			score = np.zeros( len(self.grammer.sigma) )

			for mot in self.grammer.sigma:

				t = self.normalize(mot, self.vocab)
    			#on calcule la cos similarité ente les embeddings s'ils existent
				if t in self.vocab.keys():
					embed_polyglot = self.embeddings[ self.vocab[t] ]
					score[ self.mot2id[mot] ] = cos_similarity(embedding, embed_polyglot)


			#on génère des candidats du train avec Levensthein (içi c'est word et non w car on suppose une erreur d'écriture)
			candidats, _ = self.candidats_max(word,2)
			if len(candidats) == 0 : 
				return( self.id2mot[np.argmax(score)] )

			else : 
				for mot in candidats:
					bigram_pro = self.proba_bigram(prev, mot, nex)
					unigram_pro = self.unigram[ self.mot2id[mot] ]
					score[ self.mot2id[mot] ] += coef * bigram_pro 

			return( self.id2mot[np.argmax(score)] )


    	#ne fait pas partie du polyglot
		else:
			print("erreur de grammaire...")
    		#on génère des candidats selon la distance de Levensthein
			candidats, _ = self.candidats_max(word,2)
			C = True 

    		#s'il n'y a pas de bons candidats, on genère encore plus de candidats
			if len(candidats) == 0:
				candidats, distance = self.candidats_max(word, 8)
				C = False

			#s'il n'y a toujours aucun candidat, on prend le lexique
			if len(candidats) == 0:
				candidats = self.grammer.sigma
				score = [ self.proba_bigram(prev, mot, nex)/self.unigram[self.mot2id[mot]] for mot in candidats ] 
				return( candidats[np.argmax(score)] )

			#sinon
			else:
				score = []

				#si jamais la distance max est inférieure à 2, on ne pondère pas la distance
				if C == True:
					score = [ self.proba_bigram(prev, mot, nex)/self.unigram[self.mot2id[mot]] for mot in candidats ] 

					return( candidats[np.argmax(score)] )

				#sinon on pondère
				else:
					for dist, mot in zip(distance, candidats):
						score.append( (self.proba_bigram(prev, mot, nex)/self.unigram[self.mot2id[mot]]) * (9-dist)/6 )

					return( candidats[np.argmax(score)] )























