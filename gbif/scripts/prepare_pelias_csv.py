#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for combining and formatting data so that they can be
imported into Georg. Since Georg is currently built on top of
Pelias, this means that the data can be read by the Pelias
csv-importer, see <https://github.com/pelias/csv-importer>.
"""

import pandas as pd

from ast import literal_eval


TOKENS_COLS = snakemake.params.tokens_cols
GEORG_CONFIG = snakemake.params.georg_config
PELIAS_CONFIG = snakemake.params.pelias_config


required_columns = list(set([
    GEORG_CONFIG['locationDisplayLabel'],
    PELIAS_CONFIG['name']]))


converters = dict()
tokens_dtypes = dict()
for c in TOKENS_COLS:
    converters[c] = literal_eval
    tokens_dtypes[c + 'StrLong'] = str
    tokens_dtypes[c + 'StrShort'] = str

# Load data
gbif = pd.read_table(
    snakemake.input[0], dtype=snakemake.params.gbif_dtypes)
tokens_frame = pd.read_table(
    snakemake.input[1], converters=converters, dtype=tokens_dtypes)

merged = (
    gbif.merge(tokens_frame, how='left', on='occurrenceID')
    .dropna(subset=required_columns))

# ADDENDUM JSON GBIF
addendum_json_gbif = gbif.apply(
    lambda x: x.to_json(orient='columns', force_ascii=False), axis=1)

# ADDENDUM JSON GEORG
addendum_georg = pd.DataFrame(index=merged.index)
for k, v in GEORG_CONFIG.items():
    addendum_georg[k] = merged[v]
addendum_json_georg = addendum_georg.apply(
    lambda x: x.to_json(orient='columns', force_ascii=False), axis=1)

# COMPOSE OUTPUT DATAFRAME
pelias_frame = pd.DataFrame(
    {
        'source': 'GBIF',
        'layer': snakemake.wildcards.dataset
    }, index=merged.index)

for k, v in PELIAS_CONFIG.items():
    pelias_frame[k] = merged[v]

pelias_frame['addendum_json_gbif'] = addendum_json_gbif
pelias_frame['addendum_json_georg'] = addendum_json_georg

pelias_frame.to_csv(snakemake.output[0], index=False)
