#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importación de módulos de Python necesarios para el script
import os
import pandas as pd
import numpy as np
import pylab

from Bio import SeqIO
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm #Realiza la transformació de eje a degradado de colo



#Función que realiza blast usando biopython
def makeblast(query, subj, eval, out):
	try:
		blastp_cline = ("blastp -query {0} -subject {1} -evalue {2} -outfmt "
				"'6 qseqid qseq sseqid sseq qcovs pident evalue' > {3}"
				"".format(query, subj, eval, out))
		os.system(blastp_cline)
	except:
		print('Ignorando archivo '+query+': fallo al hacer BLAST. '
		      '¿Ha comprobado el formato de input?')
		pass

#Filtro de identidad y coverage para un blast_output
def blast_filter(entrada, ident, cover, salida):
	try:
		filtrado = pd.read_csv(entrada, sep='\t', \
		names=['qseqid', 'qseq', 'sseqid', 'sseq', 'qcovs', 'pident', 'evalue'])
		filtrado = filtrado[(filtrado['qcovs'] > float(ident))
			   & (filtrado['pident'] > float(cover))]
		#Guardo el filtrado que se presenta como resultado
		filtrado.to_csv(salida, sep='\t', index = False)
		#Y finalmente detecto los rangos de coverage e identidad
		mincov = filtrado['qcovs'].min()
		maxcov = filtrado['qcovs'].max()
		minid  = filtrado['pident'].min()
		maxid  = filtrado['pident'].max()
		blasthits = len(filtrado['qcovs'])
		#Que voy a devolver para que pueda usarlos el script que llame al módulo
		return(mincov, maxcov, minid, maxid, blasthits)
	except:
		print('Ignorando archivo '+entrada+': fallo al filtrar')
		pass

#Función que calcula un histograma 2D, aka heatmap, desde un blast.tsv
def heatmap(blasttoplot, outputpng):
	try:
		(H, eje_x, eje_y) = np.histogram2d(blasttoplot.qcovs, blasttoplot.pident, bins=20)
		im = plt.imshow(H, cmap=plt.cm.Blues, norm=LogNorm(),
			extent=[eje_x[0], eje_x[-1], eje_y[0], eje_y[-1]],
			origin='lower', aspect=1)

		#Miscelánea de presentación del gráfico
		plt.title('Heatmap para los resultados de BLAST')
		plt.xlabel('Query Coverage per Subject')
		plt.ylabel('Percentage Identity')
		plt.colorbar()
		plt.savefig(outputpng)
		fig = plt.figure() #Esto me permite resetear el grafico que muestro
	except:
		print('Imposible representar el heatmap para:'+blastplot+ \
		      '. Posible fallo de BLAST')
		pass

#Función que calcula un histograma normal y corriente desde un input fasta
def histogram(fastainput, outputpng):
	try:
		tamaño = [len(record) for record in SeqIO.parse(fastainput, 'fasta')]
		pylab.hist(tamaño, bins=20)

		#Micelánea de presentación del gráfico
		pylab.title('Histograma para los resultados de BLAST')
		pylab.xlabel('Longitud de la secuencia en pb')
		pylab.ylabel('Número de secuencias')
		plt.savefig(outputpng)
		fig = plt.figure() #Esto me permite resetear el grafico que muestro
	except:
		print('Imposible representar el heatmap para:'+blastplot+ \
		      '. Posible fallo de BLAST')
		pass