# **Roteiro de Aula Prática: Análise de Bioinformática**

**Disciplina:** Biologia Molecular e Genética

**Objetivo:** Simular um fluxo de trabalho de diagnóstico molecular, desde o desenho de iniciadores (primers) até a análise comparativa de sequências gênicas, utilizando ferramentas de bioinformática online.

---

## **Etapa 1: Desenhando os Primers - O Alvo Específico**

Nesta etapa, nosso objetivo é desenhar um par de iniciadores (primers) específicos para um gene de interesse. A especificidade é crucial para garantir que nossa futura reação de PCR amplifique apenas o alvo correto.

## **Ferramenta:** [**NCBI Primer-BLAST**](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)

### **Instruções:**

1.  **Escolha UM dos três genes abaixo** para o seu estudo. Cada gene tem uma função biológica distinta.
2.  **Copie a sequência completa**, incluindo a linha que começa com `>` (cabeçalho FASTA).

#### Gene 1: Beta-globina Humana (HBB)
### *Função: Componente da hemoglobina, responsável pelo transporte de oxigênio. ###Mutações neste gene estão associadas à Anemia Falciforme.*

```
>HBB_Humano NC_000011.10:5225464-5227071
ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGC
CCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAAT
AGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGG
CTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGA
AGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCT
GCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTA
AGTTCATGTCATAGGAAGGGGATAAGTAACAGGGTACAGTTTAGAATGGGAAACAGACGAATGATTGCATCAGTGTGGAAGTCTCAGGATC
GTTTTAGTTTCTTTTATTTGCTGTTCATAACAATTGTTTTCTTTTGTTTAATTCTTGCTTTCTTTTTTTTTTCTTCTCCGCAATTTTTACT
ATTATACTTAATGCCTTAACATTGTGTATAACAAAAGGAAATATCTCTGAGATACATTAAGTAACTTAAAAAAAAACTTTACACAGTCTGC
CTAGTACATTACTATTTGGAATATATGTGTGCTTATTTGCATATTCATAATCTCCCTACTTTATTTTCTTTTATTTTTAATTGATACATAA
TCATTATACATATTTATGGGTTAAAGTGTAATGTTTTAATATGTGTACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAA
AAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAAT
GTATCATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATTTCTGCATATAAATATTTCTGC
ATATAAATTGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTTATGGTTGGGATAA
GGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGG
TCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCT
GGCCCACAAGTATCACTAA
```

#### Gene 2: Alfa-globina Humana (HBA1)
### *Função: Outro componente da hemoglobina. Mutações estão associadas às ### Talassemias Alfa.*
```
>HBA1_Humano NC_000016.10:224699-225539
TCGGCGCCCCGGCCCGGTTCTCCGCGCCGGGGCCCCGCGCCCGGCGCCCGCAAGCATAAACCCTGGCGCGCGCCGGCGCCGGGCCGGCGGC
CGCTCCGGCGCTCCGGGCCCGAGCGCCCGGCTCCCGGCCCCGGGCCCGGCGCCCGGCCTCCGCGCCGCCGTCCCACTTCCCGCCCGCCGCA
CCCGCCCCGCGCTCCCGCTCCCGCTCCCGCCCCGGCTCCCGCTCCCGCTCTCCGGCCTCCCTCCCCGCCCCCGGACCCACCCAGGGCCTTG
AGAAGC ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAG
GCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCCCACTTCGACCTGAGCCACGGCTCTGCCCAAGGTTAAGGGC
CACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACG
CGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCAC
CCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAAGCTGGAGCCTCGGTGGCC
ATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGC
```

#### Gene 3: Actina Beta Humana (ACTB)
### *Função: Proteína estrutural, componente do citoesqueleto. É um gene "housekeeping", expresso em quase todas as células.*
```
>ACTB_Humano NC_000007.14:5526938-5530432
CGGGACCTGACAGACTACCTCATGAAGATCCTCACCGAGCGCGGCTACAGCTTCACCACCACGGCCGAGCGGGAAATCGTGCGTGACATTA
AGGAGAAGCTGTGCTACGTCGCCCTGGACTTCGAGCAAGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTACGAGCTGCC
TGACGGCCAGGTCATCACCATTGGCAATGAGCGGTTCCGCTGCCCTGAGGCACTCTTCCAGCCTTCCTTCCTGGGCATGGAGTCCTGTGGC
ATCCACGAAACTACCTTCAACTCCATCATGAAGTGTGACGTGGACATCCGCAAAGACCTGTACGCCAACACAGTGCTGTCTGGCGGCACCA
CCATGTACCCTGGCATTGCCGACAGGATGCAGAAGGAGATCACTGCCCTGGCACCCAGCACAATGAAGATCAAGATCATCGCTCCCCTGAG
CGTGGCTACTCCTTCACCACCACCGCCGAGCGGGAAATCGTGCGCGAcattaaggagaagctgtgctacgtcgccctggacttcgagcaag
agatggccacggctgcttccagctcctccctggagaagagctacgagctgcctgacggccaggtcatcaccattggcaatgagcggttccg
ctgccctgaggcactcttccagccttccttcctgggcATGGAGTCCTGTGGCATCCACGAAACTACCTTCAACTCCATCATGAAGTGTGAC
GTGGACATCCGCAAAGACCTGTACGCCAACACAGTGCTGTCTGGCGGCACCACCATGTACCCTGGCATTGCCGACAGGATGCAGAAGGAGA
TCACTGCCCTGGCACCCAGCACAATGAAGATCAAGATCATCGCTCCCCTGAGCGTGGCTACTCcttcaccaccaccgccgagcgggaaatc
gtgcgcga
```

3.  Cole a sequência escolhida no campo **"PCR Template"** do site Primer-BLAST.
4.  Na seção **"Primer Parameters"**, ajuste o campo **"Product size"** para um Mínimo de `150` e um Máximo de `500`.
5.  Na seção **"Primer Pair Specificity Checking Parameters"**, certifique-se que o campo **"Organism"** está preenchido com `Homo sapiens`.
6.  Clique no botão azul **"Get Primers"** no final da página.
7.  Aguarde o resultado. A ferramenta irá sugerir vários pares de primers.

    > ⚠️ **ANOTE OS DADOS DO "PRIMER PAIR 1"**
    > * Sequência do **Forward primer**
    > * Sequência do **Reverse primer**
    > * O valor de **Tm** de cada primer.

---

## **Etapa 2: A Reação em Cadeia da Polimerase (PCR) Virtual**

Agora vamos usar os primers que você desenhou para simular a amplificação do seu gene alvo.

### **Ferramenta:** [**Virtual PCR**](https://virtual-pcr.ico2s.org/pcr/#instructions)

### **Instruções:**

1.  Copie novamente a **sequência completa do gene** que você usou na Etapa 1.
2.  No site da PCR Virtual, cole a sequência no campo **"Insert your DNA sequence here"**.
3.  Cole as sequências dos seus primers **Forward** e **Reverse** nos campos correspondentes.
4.  **Defina os parâmetros do termociclador** na interface do site:
    * **Number of cycles:** `35`
    * **Denaturation temperature:** `95` °C
    * **Denaturation time:** `30` segundos
    * **Annealing temperature:** **Use a temperatura do seu primer com o Tm mais baixo e subtraia 5°C.** (Ex: Se o menor Tm for 60°C, use 55°C).
    * **Annealing time:** `30` segundos
    * **Extension temperature:** `72` °C
    * **Extension time:** `60` segundos
5.  Clique em **"Run"** e observe o resultado no gel de agarose virtual.

**ANOTE O RESULTADO**
    
---

## **Etapa 3: Análise Comparativa e Controle de Qualidade**

Finalmente, vamos comparar as sequências dos três genes para entender o quão semelhantes eles são. Isso nos ajuda a entender a especificidade dos primers e as relações evolutivas entre os genes.

### **Ferramenta:** [**Clustal Omega**](https://www.ebi.ac.uk/jdispatcher/msa/clustalo?stype=dna)

### **Instruções:**

1.  Para cada alinhamento abaixo, copie as sequências completas (incluindo o cabeçalho `>`) dos genes indicados para a caixa de texto do Clustal Omega.
2.  Certifique-se que a opção **"DNA"** está selecionada e clique em **"Submit"**.
3.  Analise a aba **"Alignment"** nos resultados. A linha inferior mostra a conservação: um asterisco (`*`) significa que a base é idêntica em ambas as sequências naquela posição.

    ---
### **Alinhamento 1: Genes da mesma família**
* Compare: **HBB** vs. **HBA1**
* *Observe o grau de similaridade. Eles compartilham um ancestral comum?*

    ---
### **Alinhamento 2: Genes não relacionados**
* Compare: **HBB** vs. **ACTB**
* *Compare este resultado com o anterior. A similaridade é maior ou menor?*

    ---
### **Alinhamento 3: Genes não relacionados**
* Compare: **HBA1** vs. **ACTB**
* *Confirme a observação do alinhamento anterior.*

    ---

    
