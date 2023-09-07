# insertion_finder - element insertion finder in a genome through a BLAST search

Insertion_finder is a tool in Python that performs similarity searches against genomic nucleotide sequences and automatically analyzes the results, identifying occurrences without the presence of the element. By comparing the sequences with and without the element, it is possible to accurately determine the 5' and 3' insertion points. 

##   Instalation

Insertion_finder does not need to be installed. The user should only download the insertion_finder.py file.

## Requirements


## Usage
```
python insertion_finder.py -q <query file> -run 'local' -d <database file> 
python insertion_finder.py -q <query file> -run 'web'
python insertion_finder.py -q <query file> -run 'local' -d <database file> -tab <BLASTn table file>
python insertion_finder.py -q <query file> -tab <BLASTn table file> 
```
### Mandatory parameters:
```
-q <file name>      Sequence to search with (fasta or multifasta file)
-run <local|web>    Choice of running local or web BLAST search
-d <file name>      Database to BLAST against (multifasta file)
```

### Optional parameters:
```
-conf <file name>   Configuration file
-tab <file name>    BLASTn search result table (fields: qseqid,sseqid,qcovs,qlen,slen,qstart,qend)(table file)
-org <integer>      Taxid(s) to restrict the database of the BLASTn search
-out <path|name>    Output directory
-minlen <integer>   Minimum element's length in base pairs(bp) (default: 5000)
-maxlen <integer>   Maximum element's length in base pairs(bp) (default: 50000)
-mincov <integer>   Minimum % query coverage per subject (default: 30)
-maxcov <integer>   Maximum % query coverage per subject (default: 90)
-enddist <integer>  Maximum distance between block tip and query tip in base pairs(bp) (default: 50)
-cpu <integer>      Number of threads to execute the blastn search (default: 10)
-color <string>     Element RGB color that is shown by the feature table, three integers between 0 and 255 separated by commas (default: 255,0,0)
``` 

## Contact

To report bugs, to ask for help and to give any feedback, please contact Arthur Gruber (argruber@usp.br) or Giuliana L. Pola (giulianapola@usp.br).

## Versions

### 2.3.0
- Addition of the routine that checks the BLAST table and, in case some sequence of the query file is missing, performs another BLAST search with the unprocessed or incomplete sequences and repeats the routine until no sequence is missing
- Creation of the function missingquery(BLAST table file,query file in fasta format,id of the sequences that were processed) that gathers the content of the BLAST tables in a single file "blastn.tab", checks if the sequences of the query file are in the BLAST table, and if not, creates a new query file with the not found or incomplete sequences, and returns the id of the processed sequences
- Parameter run (local|web) is now mandatory
- Parameter d (database file) becomes mandatory if the run parameter is local

### 2.2.5
- The commands that perform the blast search have been placed in the blast(BLAST search parameters,query file,BLAST search output file) function to be used more than once within the code if needed
- The commands that check if the output folder exists and if it does generate a new name were put in the function rename(number,path,'dir'|'file')

### 2.2.4
- Creation of file.log right after argument validation
- Addition of screen messages displayed after argument validation in file.log

### 2.2.3
- Corrected message formatting "Invalid element, % query coverage less than valid coverage!" in file.log
- Fix error repeating cases of "no valid hits" in elements.txt and file.log

### 2.2.2
- Fix bug when getting the size of the query

### 2.2.1
- added BLAST search runtime to the files "file.log" and "elements.txt
- modification of the message when the element is invalid to make it clearer: "Invalid element, smaller than valid size!" or "Invalid element, larger than valid size!"

### 2.2.0
- added megablast option to BLASTn search (task = 'megablast')
- bug fixing, problem with maxlen and minlen parameters
- added -version parameter that shows the program version

### 2.1.1
- number of threads only appears in log and elements in local version (-run'local')
- modification of org parameter for two or more organs (single quotes were missing inside double quotes in blast's entrez_query parameter)

### 2.1.0
- o to out parameter changed
- added and validated org, taxid parameter to restrict BLASTn search database

### 2.0.1
- bug fixed (query ARBZ01000001_1 web version)

### 2.0.0
- run parameter added with local or web options
- error: choose subject twice (query ARBZ01000001_1 web version)

### 1.5.1
- fix "IndexError: list index out of range" error

### 1.5.0
- tab parameter validation
- change default maximum distance (enddist) between block and query to 50
- added and validated parameter cpu, number of processors to blast 
- changed parameter c to color ⇒ ambiguity with parameter cpu
- change in the validation of the parameters mincov and maxcov ⇒ between 0 and 100 and different from each other
- change in the validation of minlen and maxlen parameters ⇒ greater than 0 and different from each other
- addition of the duration time and parameters used in the file.log and elements.txt header

### 1.4.0
- sort hits by % coverage per subject
- changed name of min and max parameters to minlen and maxlen, minimum and maximum element size
- change default max element size (max) to 50,000
- added parameter enddist, maximum distance between the end of the block and the query so that the element is accepted when it has only 1 block
- addition of the parameters mincov and maxcov, minimum and maximum % of query coverage
- validation of the alignment for 1 block
- validation of the maximum distance between the block and the query when there is only 1 block
- creation of folder-folder-name_2 when the folder specified in parameter out already exists
- change of default minimum size (minlen) to 4000
- addition of parameter tab, table with BLASTn search results

### 1.3.1
- correction of parameter "c" validation

### 1.3.0
- fixed error message "A value is trying to be set on a copy of a slice from a DataFrame"
- fix validating query and database files using fasta format check from "SeqIO.parse" command
- fix extracting the element sequence from the query file in fasta format ⇒ search within the query file
- add hit and subject data in "file.log
- addition of the BLASTn search validation and warning in case of error
- add parameter "c", the RGB color of the element in the feature table
- error in the validation of parameter "c"

### 1.2.1
- parameters validation
- verification of the existence of the query and database files and the output directory, parameters "query", "db" and "out
- addition of the warnings in the "file.log
- creation of "output_dir" if the parameter "o" is not informed
- addition of the column "valid" in the table "elements"

### 1.2.0
- adding the data type of each parameter in the help menu
- creation of the output directory if it doesn't exist with the command "os.mkdir", when the parameter "out" is entered 
- creation of the files "blastn.tab" and "elements.txt" inside the output directory with the "open" and "os.path.join" commands
- creation of a directory for each query with the "os.mkdir" command
- creating the element's fasta and feature table files inside the folder with the query name with the commands "open" and "os.path.join
- changing the parameter names from "out", "db", "query" to "o", "d" and "q", respectively

### 1.1.0
- adding the fields "query file", "database file" and "element length" in the header of "elements.txt" (file that shows the coordinates and size of elements) to inform the parameters used in the search
- displaying the help menu when the user does not inform any parameter
- discarding the subject if the number of blocks is equal to 1
- changing the name of the BLAST output from "BLASTn_elements.txt" to "blastn.tab
