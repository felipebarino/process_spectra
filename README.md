# process_spectra

This is the repo for the process_spectra package. It is a package designed to process optical spectra of fiber optic sensors based on long period gratings (LPGs). The documentation is written in portuguese, since the project was conceived to improve research at a brazilian university lab (LITel - UFJF). If the contents of the library may be useful to you, and you do not speak portuguese, please send us an e-mail or open an issue.

Pacote python feito com o intuito de processar o espectro óptico de sensores ópticos a fibra baseados em grades de periodo longo (LPGs). Com esse pacote é possível fazer a extração de dados de conjuntos grandes de espectros, seguindo uma rotina específica. 

## Instalação:

O pacote foi colocado no PyPi, logo é possível instalar pelo pip:

```
pip install process_spectra
```

## Como usar:

Para usar o pacote, basta criar um objeto da classe *MassSpectraData*, adicionar passos com os devidos argumentos e rodar. Como um exemplo simples que extrai os vales ressonantes de espectros na pasta *spectra*:

``` python
import os
import process_spectra as ps


# Criando a lista de caminhos
files = os.listdir('spectra')
files_complete = [os.path.join('spectra', x) for x in files]

# Criando o objeto, passando a lista e o nome do arquivo para 
# salvar as informações extraídas
spectra = ps.MassSpectraData(files_complete, 'resonant_wavelengths.csv')


# Adicionando o passo de extração, com um dicionário de argumentos 
# para a função do passo
spectra.add_step(ps.funcs.find_valley, {'prominence': 5, 'ignore_errors': True})

spectra.run(ignore_errors=True)

```

Vale notar que os ignore_errors são passados para evitar que o programa encerre no caso de encontrar algum. Isso é útil quando não se tem certeza da integridade de todos os espectros, visto que se um estiver corrompido, o programa pode travar nele.

Na pasta de exemplos estão alguns scripts que foram escritos para e usados em pesquisas, e mostram mais funcionalidades do pacote.

## Documentação

A documentação não foi publicada em sites, mas foi colocada no código, então é possível encontrar ela através do comando help() do python, passando como argumento uma classe, função ou módulo, ou também lendo do código direto.

Também é possível gerar a documentação através do sphinx. Basta clonar / baixar o repositórrio, abrir um cmd na pasta sphinx e rodar o comando 

```
make html
```

A documentação será gerada em html, como um site, e pode ser encontrada no caminho *sphinx/_build/html/index.html*.