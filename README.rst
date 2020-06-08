Preparation of GBIF data for Georg
==================================

This repository contains code for processing GBIF data before importing
it into `Georg <http://github.com/naturhistoriska/georg>`_,
which is a georeferencing tool built on top of `Pelias <https://pelias.io>`_.
The data processing is carried out with the workflow management system
`Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ and a few
Python scripts.

Workflow input are source Darwin Core archives obtained from
`GBIF <https://gbif.org>`_.

For the time being, two occurrence datasets are processed:

:nhrs-nrm: GBIF-Sweden, Entomological Collections (NHRS),
		   Swedish Museum of Natural History (NRM). 
		   DOI: |nbsp| `10.15468/fpzyjx <https://doi.org/10.15468/fpzyjx>`_

:s-fbo: GBIF-Sweden, Phanerogamic Botanical Collections (S).
	    DOI: |nbsp| `10.15468/yo3mmu <https://doi.org/10.15468/yo3mmu>`_


Prerequisites
-------------

* Python 3.7
* The Python libraries `pandas <https://pandas.pydata.org>`_, 
  `spaCy <https://spacy.io>`_, and
  `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_

An easy way to get Python working on your computer is to install the
free `Anaconda distribution <http://anaconda.com/download>`_.

You can install the libraries with the following command:

.. code-block:: bash

    pip install pandas snakemake spacy


Input files
-----------

Input files should be placed at the following location:
``./data/raw/{dataset}/occurrence.txt``.


Output files
------------

After executing the workflow, you should be able to find the output
files under ``./data/processed/``.


Running the workflow
--------------------

You can run workflow by entering the following on the command-line
(adjust the number of CPU cores to fit your environment):

.. code-block:: bash

    snakemake --cores 4


The file ``./config.yaml`` determines which datasets to include, and
how the datasets are processed.

Named Entity Recognition (NER) is used to extract place names from
texts. A language model that has been trained on transcripts of
mainly Swedish labels is included in this repository.

The workflow has been executed under Python 3.7 with the following
Python packages installed:

.. code-block::

	appdirs==1.4.4
	attrs==19.3.0
	blis==0.4.1
	catalogue==1.0.0
	certifi==2020.4.5.1
	chardet==3.0.4
	ConfigArgParse==1.2.3
	cymem==2.0.3
	datrie==0.8.2
	decorator==4.4.2
	docutils==0.16
	gitdb==4.0.5
	GitPython==3.1.3
	idna==2.9
	importlib-metadata==1.6.1
	ipython-genutils==0.2.0
	jsonschema==3.2.0
	jupyter-core==4.6.3
	murmurhash==1.0.2
	nbformat==5.0.6
	numpy==1.18.5
	pandas==1.0.4
	plac==1.1.3
	preshed==3.0.2
	psutil==5.7.0
	pyrsistent==0.16.0
	python-dateutil==2.8.1
	pytz==2020.1
	PyYAML==5.3.1
	ratelimiter==1.2.0.post0
	requests==2.23.0
	six==1.15.0
	smmap==3.0.4
	snakemake==5.19.2
	spacy==2.2.4
	srsly==1.0.2
	thinc==7.4.0
	toposort==1.5
	tqdm==4.46.1
	traitlets==4.3.3
	urllib3==1.25.9
	wasabi==0.6.0
	wrapt==1.12.1
	zipp==3.1.0


License
-------

The code in this repository is distributed under the
`MIT license <https://opensource.org/licenses/MIT>`_.


Author
------

Markus Englund


.. |nbsp| unicode:: 0xA0 
   :trim:
