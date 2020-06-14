#!/usr/bin/env python3

#Importación de módulos de Python necesarios para el script
import os
import pandas as pd

from Bio import Phylo
from contextlib import redirect_stdout
from Bio.Align.Applications import MuscleCommandline



#Función que alinea secuencias con MUSCLE (previo a alineamiento)
def muscle_align(input, output):
	in_file = r'{0}'.format(input)
	out_file = r'{0}'.format(output)
	muscle_cline = MuscleCommandline(input=in_file, out=out_file)
	stdout, stderr = muscle_cline()

def muscle_maketree(alignd_sec, map):
	#Como se puede ver, muy a mi pesar he usado os en vez de biopython, como acostumbro
	#Esto es así porque el paquete MuscleCommandLine está fatal documentado,
	#y el maketree no sale por ningún lado
	z=('muscle -maketree -in {0}  -out {1}  -cluster '
	   'neighborjoining >/dev/null 2>&1  '.format(alignd_sec, map))
	os.system(z)

#Función que representa un árbol usando el módulo Phylo de biopython
def Phylo_maketree(map, output_name):
	with open(output_name, 'w') as arbol:
		with redirect_stdout(arbol):
			tree = Phylo.read(map, 'newick')
			Phylo.draw_ascii(tree)
