#!/usr/bin/env python3

#Importación de módulos de Python necesarios para el script
import os
import pandas as pd
import numpy as np
import pylab

from Bio import SeqIO
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm #Realiza la transformació de eje a degradado de colo

from Bio.Blast.Applications import NcbiblastpCommandline 



#Función que genera la base de datos
def makeblastdb(multifasta):
	creadb = 'makeblastdb -in {0} -parse_seqids -dbtype prot \
	-out basedatos >/dev/null 2>&1'.format(multifasta)
	os.system(creadb)

#Función que realiza blast usando biopython
def makeblast(multif):
	blastn_cline = NcbiblastpCommandline(query = multif, db = 'basedatos', evalue = 0.0001,
	outfmt = '6 qseqid qcovs pident evalue sseqid sseq', out = 'resultado_blast.tsv') 
	stdout, stderr = blastn_cline()

#Filtro de identidad y coverage para un blast_output
def blast_filter(ident, cover):
	filtrado = pd.read_csv('./resultado_blast.tsv', sep='\t', \
	           names=['qseqid', 'qcovs', 'pident', 'evalue', 'sseqid',  'sseq'])
	filtrado = filtrado[(filtrado['qcovs'] > float(ident)) \
	& (filtrado['pident'] > float(cover))]
	#Guardo el filtrado que se presenta como resultado
	filtrado.to_csv(r'./resultado_blast.tsv', sep='\t', index = False)
	#Y finalmente detecto los rangos de coverage e identidad
	mincov = filtrado['qcovs'].min()
	maxcov = filtrado['qcovs'].max()
	minid  = filtrado['pident'].min()
	maxid  = filtrado['pident'].max()
	#Que voy a devolver para que pueda usarlos el script que llame al módulo
	return(mincov, maxcov, minid, maxid)

#Función que calcula un histograma 2D, aka heatmap, desde un blast.tsv
def heatmap(blasttsv):
        (H, eje_x, eje_y) = np.histogram2d(blasttsv.qcovs, blasttsv.pident, bins=20)
        im = plt.imshow(H, cmap=plt.cm.Blues, norm=LogNorm(),
                        extent=[eje_x[0], eje_x[-1], eje_y[0], eje_y[-1]],
                        origin='lower', aspect=1)

        #Miscelánea de presentación del gráfico
        plt.title('Heatmap para los resultados de BLAST')
        plt.xlabel('Query Coverage per Subject')
        plt.ylabel('Percentage Identity')
        plt.savefig('heatmap.png')

#Función que calcula un histograma normal y corriente desde un input fasta
def histogram(fastainput):
        tamaño = [len(record) for record in SeqIO.parse(fastainput, 'fasta')]
        pylab.hist(tamaño, bins=20)

        #Micelánea de presentación del gráfico
        pylab.title('Histograma para los resultados de BLAST')
        pylab.xlabel('Longitud de la secuencia en pb')
        pylab.ylabel('Número de secuencias')
        plt.savefig('histograma.png')
