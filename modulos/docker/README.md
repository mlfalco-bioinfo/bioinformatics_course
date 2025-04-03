# **Uso de Containers com Docker**  

O **Docker** é uma tecnologia essencial para a bioinformática moderna, permitindo que softwares e pipelines sejam executados de forma isolada e reprodutível, sem necessidade de instalação complexa.  


## **O que são Contêineres?**  

Um **contêiner** é um ambiente virtualizado que contém tudo o que um software precisa para rodar, incluindo:  
1. Sistema operacional base leve.  
2. Dependências e bibliotecas necessárias.  
3. Código e configurações específicas do software.  

O Docker permite criar e gerenciar contêineres de maneira eficiente. Ele resolve problemas como:  
- Reprodutibilidade: Garante que o software funcione da mesma forma em qualquer sistema.  
- Facilidade de instalação: Elimina conflitos entre dependências de diferentes programas.  
- Portabilidade: Um contêiner pode ser executado em qualquer máquina com Docker instalado.  

Na bioinformática, muitos programas são distribuídos via contêineres, como **FastQC, Samtools, BWA, Kraken2, Nextflow** e **Snakemake**.  

## **Instalação do Docker**  

###  No Linux (Debian/Ubuntu)  

```bash
sudo apt update
sudo apt install -y docker.io
```

Adicione seu usuário ao grupo `docker` para rodar comandos sem `sudo`:  

```bash
sudo usermod -aG docker $USER
newgrp docker
```

Verifique a instalação:  

```bash
docker --version
```

###  No Fedora  

```bash
sudo dnf install -y docker
sudo systemctl enable --now docker
```

###  No macOS e Windows  

Baixe e instale o **Docker Desktop**:  
[![Docker Desktop](https://img.shields.io/badge/Docker-Desktop-blue?logo=docker)](https://www.docker.com/products/docker-desktop/)



##  **Comandos Básicos do Docker**  

###  Baixar e rodar um contêiner  

```bash
docker run hello-world
```

Isso baixa e executa um contêiner de teste, verificando se a instalação está funcionando.  

###  Listar imagens locais  

```bash
docker images
```

###  Baixar uma imagem sem rodar  

```bash
docker pull debian:bookworm-slim
```

###  Rodar um contêiner interativo  

```bash
docker run -it ubuntu bash
```

Isso inicia um contêiner do Ubuntu com um terminal interativo.  

###  Parar um contêiner  

```bash
docker stop <ID_do_contêiner>
```

###  Remover um contêiner ou imagem  

```bash
docker rm <ID_do_contêiner>
docker rmi <ID_da_imagem>
```

##  **Execução de Pipelines Bioinformáticos com Docker**  

Docker é amplamente usado para rodar ferramentas bioinformáticas. Veja exemplos práticos:  

###  **Rodando FastQC com Docker**  

```bash
docker run --rm -v $(pwd):/data biocontainers/fastqc:v0.11.9_cv8 fastqc /data/meu_arquivo.fastq
```

**Explicação:**  
- `--rm`: Remove o contêiner após a execução.  
- `-v $(pwd):/data`: Monta a pasta atual dentro do contêiner.  
- `biocontainers/fastqc:v0.11.9_cv8`: Imagem do FastQC no Docker Hub.  
- `fastqc /data/meu_arquivo.fastq`: Comando a ser executado dentro do contêiner.  

###  **Rodando Samtools com Docker**  

```bash
docker run --rm -v $(pwd):/data biocontainers/samtools:v1.15.1 samtools view /data/meu_arquivo.bam
```

###  **Rodando um pipeline Nextflow com Docker**  

O Nextflow pode executar workflows diretamente em contêineres:  

```bash
nextflow run nf-core/rnaseq -profile docker
```

##  **Criando seus Próprios Contêineres**  

Você pode criar um contêiner personalizado para seu pipeline usando um **Dockerfile**.  

### **Exemplo de um Dockerfile para um ambiente bioinformático**  

Crie um arquivo chamado `Dockerfile`:  

```Dockerfile
FROM debian:bookworm-slim

RUN apt update && apt install -y \
    fastqc \
    samtools \
    && rm -rf /var/lib/apt/lists/*

CMD ["bash"]
```

Construa a imagem:  

```bash
docker build -t meu_pipeline .
```

Execute o contêiner:  

```bash
docker run -it meu_pipeline
```

