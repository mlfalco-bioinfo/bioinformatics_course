# **Acesso a Bancos de Dados Públicos**  

O acesso a bancos de dados públicos é essencial na bioinformática, pois permite a obtenção de dados genômicos, transcriptômicos e metagenômicos para análises diversas. Esses bancos oferecem repositórios gratuitos de informações biológicas que podem ser utilizadas para pesquisas em evolução, medicina, biotecnologia e muitas outras áreas.  

## **Principais Bancos de Dados**  

### **NCBI (National Center for Biotechnology Information)**  
O **NCBI** é um dos maiores repositórios de dados biológicos do mundo. Ele mantém diversas bases, incluindo:  
- **GenBank**: Sequências genômicas depositadas por pesquisadores.  
- **RefSeq**: Banco de sequências de referência altamente curadas.  
- **PubMed**: Repositório de artigos científicos.  
- **Protein Data Bank (PDB)**: Estruturas tridimensionais de proteínas.  
- **BLAST (Basic Local Alignment Search Tool)**: Ferramenta de alinhamento de sequências.  

Acesse: [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)  

####  **Como baixar sequências do NCBI via terminal**  
Podemos utilizar a ferramenta **Entrez Direct (EDirect)** para interagir com o NCBI pelo terminal:  

```bash
# Instalar Entrez Direct (caso não esteja instalado)
sh -c "$(wget -qO- https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# Buscar IDs de sequências do gene BRCA1 no banco de nucleotídeos
esearch -db nucleotide -query "BRCA1[GENE] AND Homo sapiens[ORGANISM]" | efetch -format fasta > brca1_sequences.fasta
```
Isso baixará as sequências FASTA do gene **BRCA1** em humanos.


### **SRA (Sequence Read Archive)**  
O **SRA** (Sequence Read Archive) é um banco de dados mantido pelo NCBI que armazena **dados de sequenciamento de nova geração (NGS)**. Ele contém **leituras cruas** geradas por sequenciadores como Illumina, Oxford Nanopore e PacBio.  

Acesse: [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)  

####  **Baixando dados do SRA com SRA Toolkit**  
Para acessar e baixar dados do SRA, usamos a ferramenta **SRA Toolkit**:  

```bash
# Instalar SRA Toolkit
conda install -c bioconda sra-tools  

# Baixar dados de um experimento específico (exemplo: SRR000001)
prefetch SRR000001  

# Converter para FASTQ
fastq-dump --split-files SRR000001
```
Isso salvará as leituras do experimento em formato **FASTQ**, prontas para análise.


###  **Ensembl**  
O **Ensembl** é um banco de dados europeu que fornece **anotações genômicas de referência** para eucariotos. Ele inclui:  
- **Genomas completos** de diversas espécies.  
- **Variantes genéticas** e suas associações clínicas.  
- **Genes ortólogos e parálogos** entre espécies.  
- **Ferramentas de navegação e download de dados via API**.  

Acesse: [https://www.ensembl.org/](https://www.ensembl.org/)  

####  **Baixando genomas do Ensembl**  
Podemos utilizar **wget** para baixar genomas diretamente:  

```bash
# Baixar o genoma de Homo sapiens do Ensembl
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
Isso salvará o genoma de referência humano compactado (**.gz**), que pode ser utilizado para mapeamento e anotações.

##  **Ferramentas para Obtenção e Manipulação de Dados Biológicos**  

Além das ferramentas já mencionadas, existem diversas outras para manipular e processar os dados baixados. Algumas das mais importantes são:  

| **Ferramenta**  | **Descrição**  | **Comando de Instalação**  |
|----------------|---------------|-----------------------------|
| **wget**  | Baixar arquivos da web via HTTP, HTTPS ou FTP  | `sudo apt install wget` |
| **curl**  | Ferramenta alternativa para download de arquivos via terminal  | `sudo apt install curl` |
| **jq**  | Processamento de arquivos JSON (útil para APIs de bancos de dados)  | `sudo apt install jq` |
| **samtools**  | Manipulação de arquivos SAM/BAM/CRAM (sequências alinhadas)  | `conda install -c bioconda samtools` |
| **fastq-dump**  | Extração de leituras de arquivos SRA para formato FASTQ  | `conda install -c bioconda sra-tools` |

### **Exemplo de Pipeline**
Aqui está um exemplo de um fluxo de trabalho para obter, converter e analisar dados de um experimento de sequenciamento de RNA-Seq:  

```bash
# 1. Baixar os dados do SRA
prefetch SRR000001  

# 2. Converter para FASTQ
fastq-dump --split-files SRR000001  

# 3. Verificar a qualidade das leituras
fastqc SRR000001_1.fastq SRR000001_2.fastq  

# 4. Realizar o alinhamento com o genoma de referência (usando HISAT2)
hisat2 -x ref_genome -1 SRR000001_1.fastq -2 SRR000001_2.fastq -S output.sam  

# 5. Converter para formato BAM e ordenar
samtools view -bS output.sam | samtools sort -o output_sorted.bam  
```
Esse fluxo de trabalho cobre desde a obtenção dos dados até a conversão e análise.


## **Conclusão**  
O acesso a bancos de dados públicos é um dos pilares da bioinformática moderna. NCBI, SRA e Ensembl fornecem uma riqueza de informações genômicas e transcriptômicas, essenciais para pesquisas em diversas áreas da biologia. O uso eficiente de ferramentas de linha de comando para baixar e manipular esses dados é fundamental para análises bioinformáticas avançadas.  
