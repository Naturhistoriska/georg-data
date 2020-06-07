#!/usr/bin/env python
# -*- coding: utf-8 -*-

import snakemake


configfile: 'config.yaml'


rule all:
    input: expand('data/processed/{dataset}.csv', dataset=config.keys())


rule filter_sweden:
    input: 'data/raw/{dataset}/occurrence.txt'
    output: 'data/interim/filtered/{dataset}.tsv'
    params:
        dtypes = (lambda w: config[w.dataset]['columns']),
        distinct_column_set = (
            lambda w: config[w.dataset]['distinctColumnSet']),
    script: 'scripts/filter_sweden.py'


rule extract_loc_tokens:
    input: 'data/interim/filtered/{dataset}.tsv'
    output: 'data/interim/tokens/{dataset}.tsv'
    params:
        ner_dir = (
            lambda w: config[w.dataset]['tokenExtraction']['nerModelDir']),
        source_column = (
            lambda w: config[w.dataset]['tokenExtraction']['sourceColumn']),
        long_str_sep = (
            lambda w: config[w.dataset]['tokenExtraction']['longStrSep']),
        short_str_sep = (
            lambda w: config[w.dataset]['tokenExtraction']['shortStrSep'])
    script: 'scripts/extract_loc_tokens.py'


rule prepare_pelias_csv:
    input:
        'data/interim/filtered/{dataset}.tsv',
        'data/interim/tokens/{dataset}.tsv'
    output:
        'data/processed/{dataset}.csv'
    params:
        dwc_dtypes = (lambda w: config[w.dataset]['columns']),
        name = (lambda w: config[w.dataset]['peliasName']),
        displayLabel = (lambda w: config[w.dataset]['locationDisplayLabel'])
    script: 'scripts/prepare_pelias_csv.py'
