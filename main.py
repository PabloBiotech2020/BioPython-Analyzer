##!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------- ¡WARNING! --------------------------- #
# Este script requiere la instalación previa de:
# - Python3    - Blast    - Muscle

#Los siguientes módulos de Python pueden requerir instalación desde PIP:
# - Bio   -Matplotlib   - Pandas  - Termcolor 	- Progress 

#Importación de módulos de Python necesarios WWpara el script
import pandas as pd
import sys
import re
import os

from os import path
from pathlib import Path
from Bio import SeqIO
from shutil import copyfile
from termcolor import colored
from Bio import SeqIO

import blast
import muscle
import prosite

# ************************* FUNCTIONS MODULE ************************* #

# ------------ BEGGINING FUNCTIONS ------------ #

#Función de ayuda que resume el uso del script
def usage():
	print('Para usar este script, debes llamar al script usando como argumento un archivo '
	      ' multifasta query sobre el que se busca y un directorio que contenga los genbank' 
	      ' subject, además de, en este orden, los criterios de coverage, identity y e-value.'
	      ' También debe haber un archivo prosite.dat en la carpeta de ejecución ')
	sys.exit(1)

#Función que muestra la licencia
def licencia():
	print('© 2020 Pablo Ignacio Marcos López')
	print()
	print('This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.')
	print()
	print('This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.')
	print()
	print('You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.')
	print()
	print('Coded with GNU Nano')
	sys.exit(0)

def help():
	print('Recuerda, este script requiere, en el siguiente orden:')
	print('*query: Un directorio que contiene los genbank de nt para las proteínas'
	      'sobre las que se va a ejecutar la búsqueda\n*subject: El sujeto sobre el'
	      ' que se busca, que debe introducirse como (multi)fasta')

	print('*coverage: El valor mínimo (o cut-off) para la identidad, expresado como \%,'
	      'y en forma NN.NN y exclusive \n*identity: El valor mínimo (o cut-off) para el'
	      'coverage, expresado como %, y en forma NN.NN y exclusive\n*e-value: El e-value'
	      'que usará blast para filtrar los results')

	print('*Prosite.dat: Una base de datos de prosite que debe estar en la misma carpeta '
	      'en la que se ejecutan los scripts')
	sys.exit(0)

def isthisvaluevalid(value):
	if int(value)<=0 or int(value)>=100:
		print('Error: Los valores de coverage e identity deben estar entre 0 y 100')
		sys.exit(1)

#Control de argumentos
def controldeargumentos():
	if len(sys.argv) == 1:
		print('ERROR: ', end='')
		usage()
	elif (sys.argv[1] == '-v') or (sys.argv[1] == '--verbose'):
		print('Modo verbose activado')
		if len(sys.argv) == 7:
			#Renombro las variables a nombres más llevables
			isthisvaluevalid(sys.argv[3]); isthisvaluevalid(sys.argv[4])
			verbose = True; query = sys.argv[2]; subject=sys.argv[3]
			coverage=sys.argv[4]; identity=sys.argv[5]; evalue = sys.argv[6]
			return(verbose, query, subject, coverage, identity, evalue)
		elif len(sys.argv) == 4:
			print('Error: No se han proporcionado los argumentos necesarios')
			print('Seleccionando valores por defecto: cov = 50, id = 30, eval = 0.00001')
			verbose = True; query = sys.argv[2]; subject=sys.argv[3]
			coverage=50; identity=30; evalue = 0.00001
			return(verbose, query, subject, coverage, identity, evalue)
		else:
			print('ERROR: El modo verbose requiere argumentos')
			usage()
	elif (sys.argv[1] == '-h' ) or (sys.argv[1] == '--help'):
		print('Vaya, necesitas una ayudita?')
		help()
	elif (sys.argv[1] == '--options'):
		print('Al llamar al script, tienes las siguientes opciones')
		print('-h | --help:\t muestra una ayuda')
		print('-v | --verbose:\t activa el modo verbose, pero debes añadir parámetros')
		print('-license:\t muestra la licencia')
		print('--options:\t muestra las opciones')
		usage()
	elif (sys.argv[1] == '--license') or (sys.argv[1] == '--copyright'):
        	licencia()
	elif(sys.argv[1].startswith('-')) or (sys.argv[1].startswith('--')):
		print('Error:', sys.argv[1], 'no es una opción de este programa')
		usage()
	elif len(sys.argv) == 6 :
		isthisvaluevalid(sys.argv[3]); isthisvaluevalid(sys.argv[4])
		verbose = False; query = sys.argv[1]; subject=sys.argv[2]
		coverage=sys.argv[3]; identity=sys.argv[4]; evalue = sys.argv[5]
		return(verbose, query, subject, coverage, identity, evalue)
	elif len(sys.argv) == 3:
		print('Error: No se han proporcionado los argumentos necesarios')
		print('Seleccionando valores por defectp: cov = 50, id = 30, eval = 0.00001')
		verbose = False; query = sys.argv[1]; subject=sys.argv[2]
		coverage=50; identity=30; evalue = 0.00001
		return(verbose, query, subject, coverage, identity, evalue)
	else:
		print('ERROR: ', end='')
		usage()

# ------------ PRE-PARSING FUNCTIONS ------------ #

#Función que rompe el multifasta query en muchos fastas
def breakmyfasta(multifile):
	querynames = []; subindex = 0
	multi = open (multifile, 'r')
	for record in SeqIO.parse(multi, "fasta"):
		subindex += 1
		single_n = open('./tmp/fastaquery_{}.fa'.format(subindex), 'w')
		SeqIO.write(record, single_n, 'fasta')
		querynames.append(record.id)
		single_n.close
	multi.close()
	return(subindex, querynames)

def extrae_seq_translated(file):
	subject = os.path.splitext(str(file))[0]+'.fa'
	entrada = open(file, 'r'); out = open('./tmp/subjectfasta.fa', 'a')
	for seq in SeqIO.parse(entrada, "genbank") :
		for seq_feat in seq.features:
			if seq_feat.type == 'CDS':
				try:
					out.write('>{0}@{1} \n{2}\n'.format(
						  seq_feat.qualifiers['locus_tag'][0],
						  seq.name,
						  seq_feat.qualifiers['translation'][0]))
				except: #Evita fallos por noncoding genes
					pass
	entrada.close(); out.close()

#Función que comprueba si el archivo de entrada existe; si es fasta,
#lo mantengo, y si es genebank lo transformo a fasta
def check_input(file):
	if os.path.exists(file):
		archivo_starts = open(file).read()
		if archivo_starts.startswith('>') == True:
			return(file)
		else:
			try:
				subject = extrae_seq_translated(file)
				return(subject)
			except:
				print('Error: El archivo '+file+' no tiene el formato adecuado')
				usage()
	else:
		print('¡Pero esto que es! ¡Este archivo no existe! Que lío, mejor me voy')
		sys.exit(1)

# ------------ WORK PROCESSING FUNCTIONS ------------ #

def blast_to_muscle(blastinput, nowitsmuscle):
	blasttsv_i = pd.read_csv(blastinput, sep='\t', \
		     names=['qseqid', 'qseq', 'sseqid', 'sseq', 'qcovs', 'pident', 'evalue'])

	#Y necesito meterlo en un archivo
	with open(nowitsmuscle, 'a') as f:
		f.write('>'+blasttsv_i.loc[1, 'qseqid']+'\n'+blasttsv_i.loc[1, 'qseq']+'\n')
		for k in range(1, len(blasttsv_i['sseqid'])):
			f.write('>'+str(blasttsv_i.loc[k, 'sseqid'])+'\n'+str(blasttsv_i.loc[k, 'sseq'])+'\n')
	f.close()

def graph_blast(subindex):
	print('Nota: Para que los gráficos maximicen su utilidad, sólo se recomienda '
	      'graficar secuencias grandes. Desea graficar los resultados de blast? [s/N]: ', end = '')
	decission = input()
	if (( decission == 's') or ( decission == 'Sí') or
	    ( decission == 'Si')  or ( decission == 'S')):
		guardagraficos = True
		i = 1
		while i < subindex+1:
			blasttoplot = pd.read_csv('./tmp/resultado_blast_{}.tsv'.format(i), sep='\t')
			blast.heatmap(blasttoplot, './tmp/heatmap_{}.png'.format(i))
			#La función requiere los BLAST results en fasta (es decir, muscleinput)
			blast.histogram('./tmp/muscleinput_{}.fasta'.format(i), './tmp/histograma_{}.png'.format(i))
			i = i+1
		with open('./tmp/explain.txt', 'w') as tellme:
			tellme.write('* Histograma.png muestra el número de fragmentos que da '
				      'BLAST para cada longitud de secuencia. Permite ver si, por '
				      'ejemplo, tienen una distribución dada que indique distintos '
				      'tipos de alineamiento\n\n'
				      '*Heatmap.png muestra la distribución de las secuencias según '
				      'su coverage y su identity; lo ideal sería que la mayoría '
				      'estén en máximos valores de ambas; los que estén muy alejados'
				      ' de la esquina superior dcha indican baja calidad de match')
		print('['+colored('CORRECTO', 'green')+']: Gráficos generados')
		tellme.close()
		return(guardagraficos)
	else:
		guardagraficos = False
		return(guardagraficos)

# ------------ ENDING FUNCTIONS ------------ #

def print_output(i, querynames, minid, mincov, maxid, maxcov, numerodominios, blasthits):
	print('')
	print('*******************  Resultados  *******************')
	print('Proteína Query:',querynames[i-1])
	print('****BLAST****')
	print('Rango de Coverage:\t',str(mincov[i]),'-',str(maxcov[i]))
	print('Rango de Identidad:\t',str(minid[i]),'-',str(maxid[i]))
	print('Numero de Hits:',str(blasthits[i]))
	print('****Muscle****')
	print('Árbol filogenético: Consultar carpeta de resultados')
	print('****Dominios Conservados****')
	print('Se han detectado', numerodominios[i] ,'dominios conservados')
	print('****************************************************')

def makesavedir():
	#Se solicita el directorio de guardado
	print('Por favor, especifique un nombre para el directorio de guardado.'
	      ' Este estará dentro de la carpeta de trabajo: ', end='')
	workdir = Path(input())
	#Y se chequea que no exista, claro
	if path.isdir(workdir):
		print('¡Pero bueno! ¡Si este directorio ya existe! Yo no te voy a borrar nada que '
		      'para eso ya está Windows 10')
		sys.exit(1)
	else:
		os.makedirs('./{}'.format(workdir))
		return(workdir)

#Housekeeping: muevo el output e input a workdir y borro lo que no uso
def housekeeping(subindex, guardagraficos, workdir, querynames, query):
	i = 1
	while i < subindex+1:
		#Primero, creo los directorios que voy a usar
		if subindex == 1:
			rutainput = './'+str(workdir)+'/input'
			rutaoutput = './'+str(workdir)+'/output'
			os.makedirs(rutainput); os.mkdir(rutaoutput)
		else:
			expfolder = './{0}/Query_{1}'.format(workdir, querynames[i-1])
			rutainput= './{0}/Query_{1}/input'.format(workdir, querynames[i-1])
			rutaoutput= './{0}/Query_{1}/output'.format(workdir, querynames[i-1])
			os.makedirs(expfolder); os.mkdir(rutainput); os.mkdir(rutaoutput)


		#Ahora, muevo el input a su nuevo hogar
		copyfile('{}'.format(query), '{}/query.fasta'.format(rutainput))
		copyfile('./tmp/subjectfasta.fa', '{}/subject_merged.fasta'.format(rutainput))


		#A continuación, muevo el output
		os.replace('./tmp/arbol_{}.txt'.format(i), '{0}/arbol.txt'.format(rutaoutput))
		os.replace('./tmp/resultado_blast_{}.tsv'.format(i), \
			   '{0}/resultado_blast.tsv'.format(rutaoutput))
		os.replace('./tmp/Dominios_encontrados_{}.tsv'.format(i), \
			   '{0}/dominios_encontrados.tsv'.format(rutaoutput))
		os.replace('./tmp/mapa_{}.nw'.format(i),'{0}/mapa.nw'.format(rutaoutput))
		os.replace('./tmp/muscleoutput_{}.fasta'.format(i), \
                	   '{0}/muscle_aligned.fasta'.format(rutaoutput))

		#Entre los que puede haber, por supuesto, gráficos
		if guardagraficos is True:
			os.mkdir('{0}/graficos'.format(rutaoutput))
			os.replace('./tmp/heatmap_{}.png'.format(i), \
				   './{}/graficos/heatmap.png'.format(rutaoutput))
			os.replace('./tmp/histograma_{}.png'.format(i), \
				   './{}/graficos/histograma.png'.format(rutaoutput))
			copyfile('./tmp/explain.txt', './{}/graficos/explicación_graficas.txt'.format(rutaoutput))


		#Por último, borro lo que no quiero
		os.remove('./tmp/muscleinput_{}.fasta'.format(i))
		os.remove('./tmp/fastaquery_{}.fa'.format(i))
		i += 1 # Y a hacerlo todo de nuevo

	#Como solo hay uno de estos, no los itero
	os.remove('./tmp/dominios.tsv'); os.remove('./tmp/subjectfasta.fa')
	copyfile('./prosite.dat', '{}/prosite.dat'.format(workdir))
	if guardagraficos is True: os.remove('./tmp/explain.txt')
	os.removedirs('./tmp/')
	print('['+colored('CORRECTO', 'green')+']: Se han guardado los archivos')


def resumen(subindex, querynames, minid, mincov, maxid, maxcov, numerodominios, blasthits):
	print('Tiene a su disposición un resumen de los resultados. ¿Desea consultarlo? [s/N]: ', end = '')
	decission = input()
	if (( decission == 's') or ( decission == 'Sí') or
	    ( decission == 'Si')  or ( decission == 'S')):
		print('Presione la tecla correspondiente al query deseado de entre los siguientes:')
		i = 1
		while i < subindex+1:
			print('Query:', querynames[i-1], 'Tecla:', i)
			i += 1
		print('Para todas, escriba "All"; para ninguna, escriba "None"')
		querydeseado = input()
		if querydeseado == "All" or querydeseado == "all" or querydeseado == "a":
			i = 1
			while i < subindex+1:
				print_output(i, querynames, minid, mincov, maxid, maxcov, numerodominios, blasthits)
				i += 1
		elif querydeseado == "None" or querydeseado == "none" or querydeseado == "n":
			print('El programa ha finalizado')
			sys.exit(0)
		else:
			try:
				showquery = int(querydeseado) - 1
				print_output(showquery)
			except:
				print('Input Inválido. El programa ha finalizado.')
				sys.exit(1)
	else:
		print('El programa ha finalizado')
		sys.exit(0)



# ************************* MAIN MODULE ************************* #



def main():
	verbose, query, subject, coverage, identity, evalue = controldeargumentos()

	#Comienza el programa propiamente dicho
	print('Iniciando Biopython-Alalyzer...')
	if verbose is True: print('Comprobando que los archivos sean correctos...')

	try:
		os.mkdir('./tmp')
	except:
		print('Error: NO debe haber carpetas tmp en un directorio del sistema. Por favor,'
		      ' renombre la carpeta "./tmp" y vuelva a intentarlo'); sys.exit(1)
	
	#Comprobamos los GenBank y los mergeamos
	if os.path.isdir(subject):
		if verbose is True: print('Uniendo los genbank en un solo multifasta...')
		for filename in os.listdir(subject):
			check_input('{0}{1}'.format(subject, filename))
	else:	subject = check_input(subject)

	#Separo el multifasta en muchos fasta para hacer blast
	query = check_input(query)
	subindex, querynames = breakmyfasta(query)
	print('['+colored('CORRECTO', 'green')+']: Archivos preprocesados y comprobados')

	#Primero debemos parsear el archivo prosite.dat y generar un cómodo tsv + regex
	if verbose is True: print('Escaneando el archivo PROSITE...')
	endmodule = prosite.parse_dat('prosite.dat', './tmp/dominios.tsv')
	if verbose is True: print('Convirtiendo patrones de dominios a regex...')
	dominios = prosite.dat_domainstoregex('./tmp/dominios.tsv')

	#Inicializo las listas
	mincov, maxcov, minid = [None]*(subindex+1), [None]*(subindex+1), [None]*(subindex+1)
	maxid, blasthits, numerodominios = [None]*(subindex+1), [None]*(subindex+1), [None]*(subindex+1)

	#Ahora toca la parte que entra en bucle
	i = 1
	while i < subindex+1:
		#Función que realiza BLAST
		blast.makeblast('./tmp/fastaquery_{}.fa'.format(i), './tmp/subjectfasta.fa', evalue, './tmp/resultado_blast_{}.tsv'.format(i)) 
		print('['+colored('CORRECTO', 'green')+']: BLAST',i,'de',subindex,'se ha ejecutado con éxito')

		#Función que filtra los resultados de BLAST
		if verbose is True: print('Filtrando output de BLAST',i,'con el cov e ident solicitadas')
		mincov[i], maxcov[i], minid[i], maxid[i], blasthits[i] = blast.blast_filter('./tmp/resultado_blast_{}.tsv'.format(i), coverage, identity, './tmp/resultado_blast_{}.tsv'.format(i))

		#Función que convierte cada blasttsv en fasta
		if verbose is True: print('Convirtiendo resultados a FASTA para MUSCLE:',i,'de',subindex)
		blast_to_muscle('./tmp/resultado_blast_{}.tsv'.format(i),'./tmp/muscleinput_{}.fasta'.format(i))

		#Función que alinea usando muscle
		muscle.muscle_align('./tmp/muscleinput_{}.fasta'.format(i),'./tmp/muscleoutput_{}.fasta'.format(i)) #REalineado usando muscle
		if verbose is True: print('MUSCLE ha realineado las secuencias:',i,'de',subindex)

		#Función que hace un árbol newick con MUSCLE
		muscle.muscle_maketree('./tmp/muscleoutput_{}.fasta'.format(i), './tmp/mapa_{}.nw'.format(i))

		#Función que dibuja el árbol newick usando Phylo
		muscle.Phylo_maketree('./tmp/mapa_{}.nw'.format(i), './tmp/arbol_{}.txt'.format(i)) #Phylo representa el árbol gráficamente
		print('['+colored('CORRECTO', 'green')+']: Árbol filogenético',i,'de',subindex,'guardado')

		#Código que escanea los blasttsv buscando dominios en Prosite.dat
		if i == 1:
			print('Escaneando en busca de dominos (',i,'de',subindex,'). Esto podría tomar un buen rato...')
		else:
			print('Escaneando en busca de dominos (',i,'de',subindex,')')
		numerodominios[i] = prosite.search(dominios,'./tmp/resultado_blast_{}.tsv'.format(i), './tmp/Dominios_encontrados_{}.tsv'.format(i))

		i += 1

	#Decido si plotear BLAST y, de ser así, lo guardo
	guardagraficos = graph_blast(subindex)

	#Obtemho im directorio de proyeto
	directorio = makesavedir()
	if verbose is True: print('Directorio de proyecto creado')

	#Redirijo los archivos output a donde sea adecuado
	housekeeping(subindex, guardagraficos, directorio, querynames, query)

	#Ofrezco un cómodo resumen con datos destacables
	resumen(subindex, querynames, minid, mincov, maxid, maxcov, numerodominios, blasthits)

if __name__ == '__main__':
	main()
