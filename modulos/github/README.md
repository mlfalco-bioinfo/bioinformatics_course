# **GitHub e Controle de Versão**  

O uso do **Git** e do **GitHub** é essencial para o trabalho com bioinformática, pois permite rastrear mudanças em arquivos, colaborar de forma eficiente e manter versões organizadas de scripts e análises.  

## **O que é Controle de Versão?**  

O **controle de versão** é um sistema que registra as mudanças feitas em arquivos ao longo do tempo, permitindo que você:  
- Acompanhe modificações no código-fonte ou nos dados.  
- Volte para versões anteriores caso necessário.  
- Trabalhe de forma colaborativa sem conflitos.  

O **Git** é o sistema de controle de versão mais usado no mundo, e o **GitHub** é uma plataforma baseada em Git que facilita o armazenamento e compartilhamento de repositórios online.  

## **Conceitos Fundamentais do Git**  

| **Conceito** | **Descrição** |
|-------------|--------------|
| **Repositório (repo)** | Pasta contendo arquivos e histórico de versões. |
| **Commit** | Registro de uma alteração no repositório. |
| **Branch** | Linha independente de desenvolvimento dentro do repositório. |
| **Merge** | Combinação de alterações de diferentes branches. |
| **Clone** | Cópia de um repositório remoto para o computador local. |
| **Pull** | Atualização do repositório local com as mudanças do repositório remoto. |
| **Push** | Envio de commits locais para o repositório remoto. |
| **Fork** | Cópia de um repositório para sua conta, permitindo modificações sem afetar o original. |

##  **Instalação do Git**  

Para instalar o Git no Debian e derivados, use:  

```bash
sudo apt update && sudo apt install git -y
```

Para Fedora:  
```bash
sudo dnf install git -y
```

Para macOS:  
```bash
brew install git
```

Verifique a instalação:  
```bash
git --version
```

## **Configuração Inicial do Git**  

Antes de usar o Git, configure seu nome e e-mail:  

```bash
git config --global user.name "Seu Nome"
git config --global user.email "seuemail@example.com"
```

Isso garante que seus commits serão identificados corretamente.  


##  **Criação e Gerenciamento de Repositórios**  

###  Criando um novo repositório  

```bash
git init meu_projeto
cd meu_projeto
```

Isso cria uma pasta chamada `meu_projeto` com um repositório Git vazio.  

###  Adicionando arquivos ao repositório  

```bash
echo "# Meu Projeto" > README.md
git add README.md
git commit -m "Adicionando o README"
```

O comando `git add` coloca arquivos na "área de stage", e `git commit` salva as alterações no histórico.  

###  Conectando ao GitHub  

No GitHub, crie um repositório e copie o link HTTPS. Depois, execute:  

```bash
git remote add origin https://github.com/usuario/repositorio.git
git branch -M main
git push -u origin main
```

Isso conecta o repositório local ao GitHub e envia os arquivos para lá.  


##  **Fluxo de Trabalho Comum no Git**  

###  **Clonar um repositório remoto**  

```bash
git clone https://github.com/usuario/repositorio.git
cd repositorio
```

### **Fazer mudanças e registrar commits**  

```bash
nano arquivo.txt   # Edite o arquivo
git add arquivo.txt
git commit -m "Modificação no arquivo.txt"
```

### **Enviar mudanças para o GitHub**  

```bash
git push origin main
```

### **Atualizar o repositório local com mudanças remotas**  

```bash
git pull origin main
```

## **Utilização do GitHub Codespaces**  

O **GitHub Codespaces** é um ambiente de desenvolvimento baseado em nuvem, que permite rodar um repositório com todas as dependências diretamente no navegador.  

### Como ativar um Codespace?  

1. Acesse o repositório no GitHub.  
2. Clique em `Code` → `Codespaces`.  
3. Clique em `Create codespace on main`.  

Isso iniciará um ambiente de desenvolvimento remoto baseado em **Docker**.  
