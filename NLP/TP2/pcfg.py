import numpy as np
import nltk
import nltk.tree as Tree


class PCFG():
	'''
	- self.lexicon : dictionnaire de dictionnaires. A chaque noeud terminal, on associe un mot. A chaque mot est associé la proba
	- self.sigma : liste de tous les mots
	- self.gram : dictionnaire de dictionnaires. A chaque regle non terminale,on  associe un dictionnaire de ses enfants et leur proba correspondante
	- self.N : liste de toutes les règles
	- self.count_ords : dictionnaire des mots avec en valeurs leur occurence
	Class qui va créer notre répertoire de règles de grammaire et notre lexique associé à leur proba
	'''

	def __init__(self):

		self.lexicon = {}
		self.sigma = []
		self.gram = {}
		self.N = []
		self.count_words ={}


	def lexique(self, regle):
		'''
		Input:
		- regle : règle type nltk.production

		Met à jour le self.lexicon, self.sigma, self.N et self.count_words
		'''

		rule = regle.lhs()
		mot  = regle.rhs()[0]

		#on update N
		if not rule in self.N :
			self.N.append(rule)

		#on update sigma
		if mot not in self.sigma:
			self.sigma.append(mot)

		#on update count_words
		self.count_words[mot] = self.count_words.get(mot, 0) + 1

		#on update le lexique
		if rule in self.lexicon.keys():

			#si regle ->mot existe : 
			if mot in self.lexicon[rule].keys():
				self.lexicon[rule][mot] +=1

			#si regle -> mot existe pas
			else:
				self.lexicon[rule][mot] = 1

		#si j'ai même pas la règle
		else:
			self.lexicon[rule] = {mot : 1}

	def grammar(self, regle):
		'''
		Input : 
		- regle : règle type nltk.production

		Met à jour self.gram et self.N
		'''

		#règle de gauche
		regle_l = regle.lhs()

		#règle de droite
		regle_r = regle.rhs()  #sous forme de chomsky, donc de longueur 2

		#On update N
		if not regle_l in self.N:
			self.N.append(regle_l)

		for i in regle_r:
			if not i in self.N:
				self.N.append(i)

		#on update gram
		if regle_l in self.gram.keys():

			if regle_r in self.gram[regle_l].keys():
				self.gram[regle_l][regle_r] +=1

			else:
				self.gram[regle_l][regle_r] = 1

		else:
			self.gram[regle_l] = {regle_r : 1}


	def update(self, tree):
		'''
		tree : nltk.tree

		Met à jour les règles de grammaire et le lexique
		'''

		#on parcourt les règles de grammaire de l'arbre
		for regle in tree.productions():

			#si jamais c'est un mot de vocabulaire:
			if regle.is_lexical():
				self.lexique(regle)

			#une règle de grammaire
			else:
				self.grammar(regle)

	def proba(self, dic):
		'''
		Inputs:
		-dic : dictionnaire (lexique ou  règles de grammaires)

		calcule la proba de chaque règle et du lexique
		'''

		for dico in dic.keys():
			som = np.sum(list(dic[dico].values()))

			for i in dic[dico].keys():
				dic[dico][i] /= som


	def build_pcfg(self, filepath):
		'''
		Input:
		-filepath : chemin du fichier

		créer notre pcfg
		'''
		print("création pcfg...")
		with open(filepath, "r") as txt:
			for phrase in txt:

				t = nltk.tree.Tree.fromstring(phrase)
				t.collapse_unary(collapsePOS = True,  collapseRoot = True) #on supp toutes les règles untaires
				t.chomsky_normal_form(horzMarkov = 2)

				#on met à jour les règles et le lexique
				self.update(t)

		self.proba(self.gram)
		self.proba(self.lexicon)
		print("pcfg crée")