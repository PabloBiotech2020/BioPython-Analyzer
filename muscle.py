#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importación de módulos de Python necesarios para el script
import os
import pandas as pd
import sys

from Bio import Phylo
from Bio.Align.Applications import MuscleCommandline



#Función que alinea secuencias con MUSCLE (previo a alineamiento)
def muscle_align(input, output):
	try:
		in_file = r'{0}'.format(input)
		out_file = r'{0}'.format(output)
		muscle_cline = MuscleCommandline(input=in_file, out=out_file)
		stdout, stderr = muscle_cline()
	except:
		print('Imposible alinear el archivo '+query+':'
                      '¿Ha comprobado sus valores de coverage e identity?')
		pass

def muscle_maketree(alignd_sec, map):
	#Como se puede ver, muy a mi pesar he usado os en vez de biopython, como acostumbro
	#Esto es así porque el paquete MuscleCommandLine está fatal documentado,
	#y el maketree no sale por ningún lado
	try:
		z=('muscle -maketree -in {0}  -out {1}  -cluster '
		   'neighborjoining >/dev/null 2>&1  '.format(alignd_sec, map))
		os.system(z)
	except:
		print('Ignorando el archivo '+query+': Posible fallo de muscle')
		pass

#Función que representa un árbol usando el módulo Phylo de biopython
def Phylo_maketree(map, output_name):
	try:
		original = sys.stdout
		sys.stdout = open(output_name, 'w')
		tree = Phylo.read(map, 'newick')
		arbol = Phylo.draw_ascii(tree)
		sys.stdout = original
	except:
		print('Ignorando el archivo '+query+': Posible fallo de muscle')
		pass