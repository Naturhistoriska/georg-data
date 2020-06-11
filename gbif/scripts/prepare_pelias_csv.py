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

DISPLAY_LABEL_COL = snakemake.params.displayLabel
NAME_COL = snakemake.params.name
TOKENS_COLS = snakemake.params.tokens_cols

required_columns = list(set([DISPLAY_LABEL_COL, NAME_COL]))

converters = dict()
all_tokens_cols = list()
tokens_dtypes = dict()
for c in TOKENS_COLS:
    converters[c] = literal_eval
    all_tokens_cols.append(c)
    all_tokens_cols.append(c + 'StrLong')
    all_tokens_cols.append(c + 'StrShort')
    tokens_dtypes[c + 'StrLong'] = str
    tokens_dtypes[c + 'StrShort'] = str

# Load data
dwc_frame = pd.read_table(
    snakemake.input[0], dtype=snakemake.params.dwc_dtypes)
tokens_frame = pd.read_table(
    snakemake.input[1], converters=converters, dtype=tokens_dtypes)

frame = (
    dwc_frame.merge(tokens_frame, how='left', on='occurrenceID')
    .dropna(subset=required_columns))

# ADDENDUM JSON GEORG
addendum_georg = frame[all_tokens_cols].copy()
addendum_georg['locationDisplayLabel'] = frame[DISPLAY_LABEL_COL]

addendum_json_georg = addendum_georg.apply(
    lambda x: x.to_json(orient='columns', force_ascii=False), axis=1)

# ADDENDUM JSON GBIF
addendum_json_gbif = dwc_frame.apply(
    lambda x: x.to_json(orient='columns', force_ascii=False), axis=1)

# COMPOSE OUTPUT DATAFRAME
pelias_frame = pd.DataFrame(index=frame.index)
pelias_frame = pelias_frame.assign(
    source='GBIF',
    layer=frame.institutionCode + ':' + frame.collectionCode,
    lat=frame.decimalLatitude,
    lon=frame.decimalLongitude,
    name=frame[NAME_COL],
    addendum_json_gbif=addendum_json_gbif,
    addendum_json_georg=addendum_json_georg,
)

pelias_frame.to_csv(snakemake.output[0], index=False)
