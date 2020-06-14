#!/usr/bin/env python3

# ------------------------- WARNING! --------------------------- #
# Este script requiere la instalación previa de:
# - Python3    - Blast    -Muscle

#Importación de módulos de Python necesarios para el script
import pandas as pd
import sys
import re
import os

from os import path
from pathlib import Path
from Bio import SeqIO
from shutil import copyfile
from contextlib import redirect_stdout
from termcolor import colored
from Bio import SeqIO

import blast
import muscle
import prosite


# ------------------------- BEGINNING MODULE --------------------------- #



#Función de ayuda que resume el uso del script
def usage():
	print('Para usar este script, debes llamar al script usando como argumento un archivo'
	      ' fasta query y otro multifasta subject sobre el que se forma la db de' 
	      ' búsqueda, además de, en este orden, los criterios de coverage e identity')
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

#Control de argumentos
if len(sys.argv) == 1:
	print('ERROR: ', end='')
	usage()
elif (sys.argv[1] == '-v') or (sys.argv[1] == '--verbose'):
	if len(sys.argv) == 6:
		print('Modo verbose activado')
		verbose = True; entryquery = sys.argv[2]; entrysubject=sys.argv[3]; 
		coverage=sys.argv[4]; identity=sys.argv[5]
	else:
		print('ERROR: El modo verbose requiere argumentos')
		usage()
elif (sys.argv[1] == '-h' ) or (sys.argv[1] == '--help'):
	print('Vaya, necesitas una ayudita?')
	print('Recuerda, este script requiere, en el siguiente orden:')
	print('*query: El fichero en formato fasta que contiene las proteínas sobre las que'
	      'se va a ejecutar la búsqueda\n*subject: El sujeto sobre el que se busca, es'
	      'decir, el que usará blast para formar su database. También ha de ser fasta')
	print('*coverage: El valor mínimo (o cut-off) para la identidad, expresado como %,'
	      'y en forma NN.NN y exclusive \n*identity: El valor mínimo (o cut-off) para el'
	      'coverage, expresado como %, y en forma NN.NN y exclusive')
	sys.exit(0)
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
elif len(sys.argv) != 5:
        print('ERROR: ', end='')
        usage()
else:
	#Comienza el programa propiamente dicho
	print('Este script permite comparar una proteína o una serie de ellas contra una db,'
	      'alinear los hits o coincidencias y generar un árbol filogenético N-J' 
	      'usando MUSCLE, para finalmente mostrar info de las mismas usando PROSITE')
	#Renombro las variables a nombres más llevables
	verbose = False; emtryquery = sys.argv[1]; entrysubject=sys.argv[2]; 
	coverage=sys.argv[3]; identity=sys.argv[4]

#A continuación, comprobamos si el archivo de entrada existe; si es fasta. lo mantengo,
#y si es genebank lo transformo a fasta
def check_input(file):
	archivo = Path(file)
	if archivo.is_file() :
		archivo_starts = open(archivo).read()
		if archivo_starts.startswith('>') == True:
			return(file)
		elif (file.endswith('.gb')) or (file.endswith('.genbank')):
			subject = os.path.splitext(file)[0]+'.fa'
			SeqIO.convert(file, 'genbank', subject, 'fasta')
			return(subject)
		else:
			print('Error: El archivo'+file+'no tiene el formato adecuado')
			usage()
	else:
		print('¡Pero esto que es! ¡Este archivo no existe! Que lío, mejor me voy')
		sys.exit(1)

if verbose is True: print('Comprobando que los archivos sean correctos...')
subject = check_input(entrysubject)
query = check_input(entryquery) #Función que comprueba el formato FASTA o GenBnak



# ------------------------- BLAST MODULE --------------------------- #



#Lo primero que debemos hacer es comprobar que tanto el fichero query como el subject 
#tienen formato fasta. Para ello, uso el módulo blast que he programado

blast.makeblastdb(subject) #Función que crea una base de datos de BLAST
if verbose is True: print('Creada base de datos de BLAST')

if verbose is True: print('Corriendo BLAST...')
blast.makeblast(query) #Función que realiza BLAST
print('['+colored('CORRECTO', 'green')+']: BLAST se ha ejecutado con éxito')

if verbose is True: print('Filtrando output de BLAST con el cov e ident solicitadas')
mincov, maxcov, minid, maxid = blast.blast_filter(coverage, identity) 
#Filtro del output de BLAST, teniendo cuidado de guardar los parámetros que sé que devuelve
# esta función porque me he leido su documentación en la repo de codeberg



# ------------------------- MUSCLE MODULE --------------------------- #



#A continuación, genero con el módulo muscle el realineamiento y el árbol

#Para ello, necesito filtrar las columnas que voy a usar y meter el > de fasta
if verbose is True: print('Leyendo resultados de BLAST...')
blasttsv = pd.read_csv('./resultado_blast.tsv', sep='\t')

if verbose is True: print('Convirtiendo resultados a FASTA para MUSCLE...')
#Y necesito meterlo en un archivo
with open('muscleinput.fasta', 'a') as f:
	for i in range(len(blasttsv['qseqid'])):
		f.write('>'+blasttsv.loc[i, 'qseqid']+'\n'+blasttsv.loc[i, 'sseq']+'\n')

muscle.muscle_align('muscleinput.fasta','muscleoutput.fasta') #REalineado usando muscle
if verbose is True: print('['+colored('CORRECTO', 'green')+']: MUSCLE ha realineado las secuencias')

muscle.muscle_maketree('muscleoutput.fasta', 'mapa.nw') #Genera el árbol filogenético
if verbose is True: print('El musculoso árbol ha enraizado')

muscle.Phylo_maketree('mapa.nw', 'arbol.txt') #Phylo representa el árbol gráficamente
print('['+colored('CORRECTO', 'green')+']: Árbol filogenético guardado como arbol.txt')



# ------------------------- PROSITE MODULE --------------------------- #



#Primero debemos parsear el archivo prosite.dat y generar un cómodo tsv
if verbose is True: print('Escaneando el archivo PROSITE...')
prosite.parse_dat('prosite.dat', 'dominios.tsv')

if verbose is True: print('Convirtiendo patrones de dominios a regex...')
dominios = prosite.dat_domainstoregex('dominios.tsv')

if verbose is True: print('Escaneando en busca de dominos. Esto podría tomar un buen rato...')
numerodominios = prosite.search(dominios,'resultado_blast.tsv')
if verbose is True: print('['+colored('CORRECTO', 'green')+']: Escaneo terminado')



# ------------------------- BLAST-GRAPHING MODULE --------------------------- #



print('Nota: Para que los gráficos maximicen su utilidad, sólo se recomienda'
      'graficar secuencias grandes. Desea graficar los resultados de blast? [s/N]')
decission = input()
if (( decission == 's') or ( decission == 'Sí') or 
    ( decission == 'Si')  or ( decission == 'S')):
	guardagraficos = True
	blast.heatmap(blasttsv)
	#La función requiere los BLAST results en fasta (es decir, muscleinput)
	blast.histogram('muscleinput.fasta')
	with open('explain.txt', 'w') as explainme:
		with redirect_stdout(explainme):
			print('* Histograma.png muestra el número de fragmentos que da '
			      'BLAST para cada longitud de secuencia. Permite ver si, por '
			      'ejemplo, tienen una distribución dada que indique distintos '
			      'tipos de alineamiento\n'
			      '*Heatmap.png muestra la distribución de las secuencias según '
			      'su coverage y su identity; lo ideal sería que la mayoría '
			      'estén en máximos valores de ambas; los que estén muy alejados'
			      ' de la esquina superior dcha indican baja calidad de match')



# ------------------------- ENDING MODULE --------------------------- #



#Se solicita el directorio de guardado
print('Por favor, especifique un directorio de trabajo para el output')
workdir = Path(input())
#Y se chequea que no exista, claro
if path.isdir(workdir):
	print('¡Pero bueno! ¡Si este directorio ya existe! Yo no te voy a borrar nada que '
	      'para eso ya está Windows 10')
	sys.exit(1)
else:
	rutainput = './'+str(workdir)+'/input'; rutaoutput = './'+str(workdir)+'/output'
	os.makedirs(rutainput); os.mkdir(rutaoutput)
	print('Directorio de trabajo creado')

#Housekeeping: muevo el output e input a workdir y borro lo que no uso
os.replace('arbol.txt', './{}/output/arbol.txt'.format(workdir))
os.replace('resultado_blast.tsv', './{}/output/resultado_blast.tsv'.format(workdir))
os.replace('Dominios_encontrados.tsv', './{}/output/Dominios_encontrados.tsv'.format(workdir))

if guardagraficos is True:
	rutagraphs = './'+str(workdir)+'/graficos'
	os.mkdir(rutagraphs)
	os.replace('heatmap.png', './{}/graficos/heatmap.png'.format(workdir))
	os.replace('histograma.png', './{}/graficos/histograma.png'.format(workdir))
	os.replace('explain.txt', './{}/graficos/explicación_graficas.txt'.format(workdir))

copyfile('{}'.format(query), './{}/input/query.fasta'.format(workdir))
copyfile('{}'.format(subject), './{}/input/subject.fasta'.format(workdir))
copyfile('prosite.dat', './{}/input/prosite.dat'.format(workdir))
copyfile('prosite.doc', './{}/input/prosite.doc'.format(workdir))

os.remove('mapa.nw'); os.remove('muscleinput.fasta'); os.remove('muscleoutput.fasta')
os.remove('dominios.tsv')
for filename in Path('.').glob('basedatos*'): #Con este comando MUCHO OJO
    os.remove(filename)


print('')
print('*******************  Resultados  *******************')
print('Directorio de trabajo:\t',workdir)
print('****BLAST****')
print('Rango de Coverage:\t',str(mincov),'-',str(maxcov))
print('Rango de Identidad:\t',str(minid),'-',str(maxid))
print('Archivo con los Hits disponible en /'+str(workdir)+'/output')
print('****Muscle****')
print('Árbol filogenético: Ver '+str(workdir)+'/output/arbol.txt')
print('****Dominios Conservados****')
print('Se han detectado', numerodominios ,'dominios conservados')
print('Guardados en el directorio/'+str(workdir)+'/output')
print('****************************************************')

