# insertion_finder - element insertion finder in a genome through a BLAST search
(c) 2021. Arthur Gruber & Giuliana Pola

Usage: *insertion_finder.py -q* sequence to search with *-d* database to BLAST against

**Mandatory parameters:**

* *-q* (file) Sequence to search with

* *-d* (file) Database to BLAST against

**Optional parameters:**
  
* *-o* (path)   Output directory
  
* *-min* (int)    Minimum element's length in base pairs(bp) (default: >=5000)

* *-max* (int)    Maximum element's length in base pairs(bp) (default: <=150.000)


## 06/01/2022
- adição dos campos “query file”, “database file” e “element length” no cabeçalho do “elements.txt” (arquivo que mostra as coordenadas e o tamanho dos elementos) para informar os parâmetros utilizados na busca
- exibição do menu de ajuda quando o usuário não informa nenhum parâmetro
- descarte do subject caso o número de blocos seja igual a 1
- mudança do nome da saída de BLAST de “BLASTn_elements.txt” para “blastn.tab”

## 07/01/2022
- adição do tipo de dados de cada parâmetro no menu ajuda
- criação do diretório de saída caso ele não exista com o comando “os.mkdir”, quando o parâmetro “out” é informado 
- criação dos arquivos “blastn.tab” e “elements.txt” dentro do diretório output com os comandos “open” e “os.path.join”
- criação de um diretório para cada query com o comando “os.mkdir”
- criação dos arquivos fasta e feature table do elemento dentro da pasta com o nome da query com os comandos “open” e “os.path.join”
- alteração do nome dos parâmetros de “out”, “db”, “query” para “o”, “d” e “q”, respectivamente

## 10/01/2022
- validação do parâmetros
- verificação da existência dos arquivos de query e database e do diretório output, parâmetros “query”, “db” e “out”
- adição dos avisos no “file.log”
- criação do “output_dir” caso o parâmetro “o” não for informado
- adição da coluna “valid” na tabela “elements”

## 11/01/2022
- correção da mensagem de erro “A value is trying to be set on a copy of a slice from a DataFrame”