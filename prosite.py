#!/usr/bin/env python3

#Importación de módulos de Python necesarios para el script
import re
import pandas as pd

from Bio.ExPASy import Prosite, Prodoc
from contextlib import redirect_stdout
from progress.bar import FillingCirclesBar



#Function to parse the prosite.dat file
def parse_dat(prositedat, output_file):
	with open(output_file, 'w') as dominios:
		with redirect_stdout(dominios):
			print('name\taccession\tdescription\tpattern')
			handle = open('prosite.dat','r')
			records = Prosite.parse(handle)
			for record in records:
    				print(record.name+'\t'+record.accession+'\t'
				      +record.description+'\t'+record.pattern)

#Función que convierte las expresiones de prosite a regex
def dat_domainstoregex(entry_file):
	dominios = pd.read_csv(entry_file, sep='\t') #Lectura del archivo input
	dominios.dropna(subset=['pattern'], inplace=True) #Corta filas sin dominio
	dominios.reset_index(inplace=True) #Y resetea el index para no dar problemas
	#Convierte los pattern de prosite en regex válidos
	dominios['pattern'] = dominios.pattern.replace({'-':'','x':'.','\(':'{',
							'\)':'}'}, regex=True)
	return(dominios)

#Función de búsqueda de los prosite patterns en un documento dado
#Para mayor comodidad, la función devuelve los pattern en formato regex
#También incorpora una barra de progreso porque es de ejecución larga
def search(inputlist, protein_seqs):
	numerodominios = 0 #Inicializa el total de matches
	lineaalinea = pd.read_csv(protein_seqs, sep='\t')
	bar = FillingCirclesBar('Buscando dominios...', max =
	                         len(inputlist['pattern'])*len(lineaalinea['qseqid']))
	with open('Dominios_encontrados.tsv', 'w') as foundomains:
		with redirect_stdout(foundomains):
			print('blast hit\tname\taccession\tdescription\tpattern')
			for j in range(len(lineaalinea['qseqid'])):
				for i in range(len(inputlist['pattern'])):
					busca = inputlist.loc[i, 'pattern']
					prosearch = lineaalinea.loc[j, 'sseq']
					match = re.search(busca, prosearch, flags=re.I)
					bar.next()
					if match:
						print( lineaalinea.loc[j,'qseqid']+'\t' \
						      +inputlist.loc[i, 'name']+'\t' \
						      +inputlist.loc[i, 'accession']+'\t' \
						      +inputlist.loc[i, 'description']+'\t' \
					              +inputlist.loc[i, 'pattern'])
						numerodominios += 1
	bar.finish()
	return(numerodominios)
