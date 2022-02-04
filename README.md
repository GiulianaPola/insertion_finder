# insertion_finder - element insertion finder in a genome through a BLAST search
(c) 2021. Arthur Gruber & Giuliana Pola

Usage: *insertion_finder.py -q* sequence to search with *-d* database to BLAST against *-run* 'local'

*insertion_finder.py -q* query file *-tab* BLASTn table file

*insertion_finder.py -q* query file *-run* 'web'

**Mandatory parameters:**

* *-q* (fasta or multifasta file) - Sequence to search with

* *-d* (multifasta file) - Database to BLAST against

* *-tab* (table file) - BLASTn search result table (fields: qseqid,sseqid,qcovs,qlen,slen,qstart,qend)
* *-run* (local|web) - Choice of running local or web BLAST search

**Optional parameters:**
* *-org* (int) - Taxid(s) to restrict the database of the BLASTn search

* *-out* (path) - Output directory
  
* *-minlen* (int) - Minimum element's length in base pairs(bp) (default: 5000)

* *-maxlen* (int) - Maximum element's length in base pairs(bp) (default: 50000)

* *-mincov* (int) - Minimum % query coverage per subject (default: 30)

* *-maxcov* (int) - Maximum % query coverage per subject (default: 90)

* *-enddist* (int) - Maximum distance between block tip and query tip in base pairs(bp) (default: 50)

* *-cpu* (int) - Number of threads to execute the blastn search (default: 10)

* *-color* (int) - Element RGB color that is shown by the feature table, three integers between 0 and 255 separated by commas (default: 255,0,0)

## 06/01/2022 (1.1.0)
- adição dos campos “query file”, “database file” e “element length” no cabeçalho do “elements.txt” (arquivo que mostra as coordenadas e o tamanho dos elementos) para informar os parâmetros utilizados na busca
- exibição do menu de ajuda quando o usuário não informa nenhum parâmetro
- descarte do subject caso o número de blocos seja igual a 1
- mudança do nome da saída de BLAST de “BLASTn_elements.txt” para “blastn.tab”

## 07/01/2022 (1.2.0)
- adição do tipo de dados de cada parâmetro no menu ajuda
- criação do diretório de saída caso ele não exista com o comando “os.mkdir”, quando o parâmetro “out” é informado 
- criação dos arquivos “blastn.tab” e “elements.txt” dentro do diretório output com os comandos “open” e “os.path.join”
- criação de um diretório para cada query com o comando “os.mkdir”
- criação dos arquivos fasta e feature table do elemento dentro da pasta com o nome da query com os comandos “open” e “os.path.join”
- alteração do nome dos parâmetros de “out”, “db”, “query” para “o”, “d” e “q”, respectivamente

## 10/01/2022 (1.2.1)
- validação do parâmetros
- verificação da existência dos arquivos de query e database e do diretório output, parâmetros “query”, “db” e “out”
- adição dos avisos no “file.log”
- criação do “output_dir” caso o parâmetro “o” não for informado
- adição da coluna “valid” na tabela “elements”

## 11/01/2022 (1.3.0)
- correção da mensagem de erro “A value is trying to be set on a copy of a slice from a DataFrame”
- correção da validação dos arquivos query e database usando a verificação do formato fasta a partir do comando “SeqIO.parse”
- correção na extração da sequência do elemento a partir do arquivo query no formato fasta ⇒ busca dentro do arquivo query
- adição dos dados do hit e do subject no “file.log”
- adição da validação da busca BLASTn e aviso caso erro
- adição do parâmetro “c”, a cor RGB do elemento na feature table
- erro na validação do parâmetro “c”

## 12/01/2022 (1.3.1)
- correção da validação do parâmetro “c”

## 14/01/2022 (1.4.0)
- ordenação dos hits por % de cobertura por subject
- alteração do nome dos parâmetros min e max para minlen e maxlen, tamanho mínimo e máximo do elemento
- alteração do tamanho máximo padrão do elemento (max) para 50.000
- adição do parâmetro enddist, distância máxima entre a extremidade do bloco e da query para que o elemento seja aceito quando tiver apenas 1 bloco
- adição dos parâmetros mincov e maxcov, % mínima e máxima de cobertura da query
- validação do alinhamento para 1 bloco
- validação da distância máxima entre o bloco e a query quando há só 1 bloco
- criação da pasta nome-de-pasta_2 quando a pasta informada no parâmetro out já existir
- alteração do tamanho mínimo padrão (minlen) para 4000
- adição do parâmetro tab, tabela com resultado da busca BLASTn

## 15/01/2022 (1.5.0)
- validação do parâmetro tab
- alteração da distância máxima (enddist) padrão entre o bloco e a query para 50
- adição e validação do parâmetro cpu, número de processadores para executar o blast 
- alteração do parâmetro c para color ⇒ ambiguidade com o parâmetro cpu
- alteração na validação dos parâmetros mincov e maxcov ⇒ entre 0 e 100 e diferentes entre si
- alteração na validação dos parâmetros minlen e maxlen ⇒ maiores que 0 e diferentes entre si
- adição do tempo de duração e parâmetros usados no cabeçalho do file.log e elements.txt

## 21/01/2022 (1.5.1)
- correção do erro "IndexError: list index out of range"

## 26/01/2022 (2.0.0)
- adição do parâmetro run com as opções local ou web
- erro: escolha de subject duas vezes (query ARBZ01000001_1 versão web)

## 29/01/2022 (2.0.1)
- erro corrigido (query ARBZ01000001_1 versão web)

## 30/01/2022 (2.1.0)
- mudança do parâmetro o to out
- adição e validação do parâmetro org, taxid para restringir o database da busca BLASTn

## 31/01/2022 (2.1.1)
- number of threads só aparece no log e no elements na versão local (-run’local’)
- modificação do parâmetro org para dois ou mais organismos (faltava aspas simples dentro das aspas duplas no parâmetro entrez_query do blast)

## 01/02/2022 (2.2.0)
- adição da opção megablast na busca BLASTn (task = ’megablast’)
- correção de erro, problema nos parâmetros maxlen e minlen
- adição do parâmetro -version que mostra qual a versão do programa

## 03/02/2022 (2.2.1)
- adição do tempo de execução da busca BLAST no arquivos “file.log” e “elements.txt”
- modificação da mensagem quando o elemento é inválido para ficar mais claro: “Invalid element, smaller than valid size!” ou “Invalid element, larger than valid size!"