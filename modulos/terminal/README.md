   - **Comandos básicos para navegação e manipulação de arquivos**
     - Listar arquivos e diretórios:
       ```bash
       ls -l    # Lista arquivos com detalhes
       ls -a    # Exibe arquivos ocultos
       ```
     - Navegação entre diretórios:
       ```bash
       cd diretório    # Entra no diretório
       cd ..           # Volta um nível
       ```
     - Criar e remover arquivos/diretórios:
       ```bash
       mkdir nova_pasta   # Cria um diretório
       touch arquivo.txt  # Cria um arquivo vazio
       rm arquivo.txt     # Remove um arquivo
       rmdir pasta_vazia  # Remove um diretório vazio
       rm -r pasta        # Remove um diretório e seu conteúdo
       ```
     - Mover e copiar arquivos:
       ```bash
       mv arquivo.txt novo_local/   # Move arquivo
       cp arquivo.txt copia.txt     # Copia arquivo
       ```
     - Exibir conteúdo de arquivos:
       ```bash
       cat arquivo.txt   # Exibe o conteúdo
       less arquivo.txt  # Permite rolar pelo conteúdo
       head -n 10 arquivo.txt  # Mostra as 10 primeiras linhas
       tail -n 10 arquivo.txt  # Mostra as 10 últimas linhas
       ```
     - Procurar dentro de arquivos:
       ```bash
       grep "palavra" arquivo.txt   # Busca uma palavra no arquivo
       grep -r "palavra" pasta/     # Busca recursivamente em diretórios
       ```
   
   - **Introdução ao Bash scripting**
     - Criar um script simples:
       ```bash
       echo "#!/bin/bash" > script.sh
       echo "echo Hello, World!" >> script.sh
       chmod +x script.sh  # Dá permissão de execução
       ./script.sh         # Executa o script
       ```
     - Variáveis em Bash:
       ```bash
       nome="Mateus"
       echo "Meu nome é $nome"
       ```
     - Estruturas de controle:
       ```bash
       if [ -f "arquivo.txt" ]; then
           echo "O arquivo existe!"
       else
           echo "O arquivo não existe!"
       fi
       ```
     - Loops:
       ```bash
       for i in {1..5}; do
           echo "Iteração $i"
       done
       ```
   
   - **Gerenciamento de permissões e processos**
     - Modificar permissões de arquivos:
       ```bash
       chmod 755 script.sh  # Define permissões de execução
       chmod +x script.sh   # Dá permissão de execução
       ```
     - Modificar proprietário e grupo:
       ```bash
       chown usuario:grupo arquivo.txt  # Muda dono do arquivo
       ```
     - Ver processos em execução:
       ```bash
       ps aux   # Lista processos
       top      # Exibe processos em tempo real
       htop     # Interface interativa para monitorar processos
       ```
     - Finalizar processos:
       ```bash
       kill PID  # Encerra um processo pelo ID
       pkill nome_processo  # Encerra pelo nome
       ```
