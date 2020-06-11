# BioPython-Analyzer
<img src="https://upload.wikimedia.org/wikipedia/commons/1/13/Biopython_logo.png" width="30%">
<br>
Este script permite comparar una proteína o una serie de ellas contra una db, alinear los hits o coincidencias con BLAST y generar un árbol filogenético N-J usando MUSCLE, para finalmente mostrar info de los dominios presentes en las mismas usando PROSITE

## Requisitos previos
Para ejecutar este script, es necesario tener instalado:

| Paquete | Descripción | Orden para instalar |
| -------- | -------- | -------- |
| Python 3.8 | Lenguaje de programación usado | `sudo apt install python3.8` |
| Biopython | Herramientas para Biología Computacional | `pip install biopython` |
| BLAST | Utilidad de alineamiento de secuencias | `sudo apt install ncbi-blast+` |
| MUSCLE 3.8 | Programa de alineamiento múltiple | `sudo apt install muscle` |

Alternativamente, como one-liner: `sudo apt install python3.8 ncbi-blast+ muscle; pip install biopython`

## Descarga
Para usar este método, es necesario tener instalado git en el sistema (`sudo apt install git-all`)
1. Clona este repositorio: `git clone https://codeberg.org/FlyingFlamingo/BioPython-Analyzer`
2. Entra en el directorio: `cd BioPython-Analyzer`
3. Ejecútalo

## Uso
El script consta de un módulo principal, que debe ejecutarse como 

`python3 main.py query subject coverage identity`

donde:

* `query`: El fichero **en formato fasta** que contiene las proteínas sobre las que se va a ejecutar la búsqueda
* `subject`: El sujeto sobre el que se busca, es decir, el que usará blast para formar su database. **También ha de ser fasta**
* `coverage`: El valor mínimo (o cut-off) para la identidad, expresado como %, y en forma NN.NN [y exclusive](https://dle.rae.es/exclusive)
* `identity`: El valor mínimo (o cut-off) para el coverage, expresado como %, y en forma NN.NN [y exclusive](https://dle.rae.es/exclusive)

Durante la ejecución, main.py llama a los módulos secundarios, blast.py, muscle.py y prosite.py. 

## Opciones
Al llamar al script principal, tenemos las siguientes opciones:

| Opción | Descripción |
| -------- | -------- |
| ``--help | -h`` | Muestra la ayuda |









