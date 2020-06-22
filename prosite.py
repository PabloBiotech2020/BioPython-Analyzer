#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importación de módulos de Python necesarios para el script
import re
import pandas as pd

from Bio.ExPASy import Prosite
from progress.bar import FillingCirclesBar



#Function to parse the prosite.dat file
def parse_dat(prositedat, output_file):
	try:
		with open(output_file, 'w') as dominios:
			dominios.write('name\taccession\tdescription\tpattern\n')
			handle = open(prositedat,'r')
			records = Prosite.parse(handle)
			for record in records:
    				dominios.write(	record.name+'\t'+record.accession+'\t'
					       +record.description+'\t'+record.pattern+'\n')
			handle.close()
		dominios.close()
		return(dominios)
	except:
		print('No se ha podido leer el archivo: '+prositedat+'. Abortando módulo...')

#Función que convierte las expresiones de prosite a regex
def dat_domainstoregex(entry_file):
	try:
		dominios = pd.read_csv(entry_file, sep='\t') #Lectura del archivo input
		dominios.dropna(subset=['pattern'], inplace=True) #Corta filas sin dominio
		dominios.reset_index(inplace=True) #Y resetea el index para no dar problemas
		#Convierte los pattern de prosite en regex válidos
		dominios['pattern'] = dominios.pattern.replace({'-':'','x':'.','\{':'^[','\}':']',
								'\(':'{','\)':'}'}, regex=True)
		return(dominios)
	except:
		print('No se ha podido leer el archivo: '+entry_file+'. Abortando módulo...')

#Función de búsqueda de los prosite patterns en un documento dado
#Para mayor comodidad, la función devuelve los pattern en formato regex
#También incorpora una barra de progreso porque es de ejecución larga
def search(inputlist, protein_seqs, tsvsalida):
	try:
		numerodominios = 0 #Inicializa el total de matches
		lineaalinea = pd.read_csv(protein_seqs, sep='\t')
		bar = FillingCirclesBar('Buscando dominios...', max =
	 				 len(inputlist['pattern'])*len(lineaalinea['qseqid']))
		with open(tsvsalida, 'w') as found:
			found.write('blast hit\tname\taccession\tdescription\tpattern\n')
			for j in range(len(lineaalinea['qseqid'])):
				for k in range(len(inputlist['pattern'])):
					busca = inputlist.loc[k, 'pattern']
					prosearch = lineaalinea.loc[j, 'sseq']
					match = re.search(busca, prosearch, flags=re.I)
					bar.next()
					if match:
						found.write( lineaalinea.loc[j,'qseqid']+'\t' \
						            +inputlist.loc[k, 'name']+'\t' \
						            +inputlist.loc[k, 'accession']+'\t' \
						            +inputlist.loc[k, 'description']+'\t' \
					                    +inputlist.loc[k, 'pattern']+'\n')
						numerodominios += 1
		found.close(); bar.finish()
		return(numerodominios)
	except:
		print('Fallo al buscar dominios')
		pass
