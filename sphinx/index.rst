.. process_spectra documentation master file, created by
   sphinx-quickstart on Tue Apr 13 14:18:05 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Bem vindo(a) à documentação do process_spectra!
===============================================
A biblioteca foi desenvolvida com o intuito de processar grandes quantidades de
espectros com um conjunto determinado de funções. Por padrão, ela assume o
formato de saída do optisystem. É modular e permite a criação e integração de
novas funções, desde que recebam como parâmetro o espectro e o dicionário de
informações, e retorne o novo espectro e o dicionário com as informações novas.

A função pode ter outros parâmetros, mas o espectro deve ser o primeiro, o
dicionário de informações o segundo, e outros deverão ser passados através do
dicionário 'kwargs'.

Getting started
===============
Para utilizar a biblioteca, vale a pena observar o exemplo mais elementar: A
extração do vale ressonante, em posição e atenuação. O exemplo está no github,
na pasta exemplos. Sob o nome 'get\_resonant\_info.py'.


Na pasta exemplos tem mais alguns scripts já prontos que fazem o uso das
funcionalidades básicas.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   docs


Índices e tabelas
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
