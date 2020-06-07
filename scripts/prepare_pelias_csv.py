#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

from ast import literal_eval

DISPLAY_LABEL_COL = snakemake.params.displayLabel
NAME_COL = snakemake.params.name

required_columns = list(set([DISPLAY_LABEL_COL, NAME_COL]))


# Load data
dwc_frame = pd.read_table(
    snakemake.input[0], dtype=snakemake.params.dwc_dtypes)
tokens_frame = pd.read_table(
    snakemake.input[1], converters={'locTokens': literal_eval})

frame = dwc_frame.merge(
    tokens_frame, how='left', on='catalogNumber'
).dropna(subset=required_columns)

# ADDENDUM JSON GEORG
addendum_georg = pd.DataFrame(index=frame.index)
addendum_georg['locationDisplayLabel'] = frame[DISPLAY_LABEL_COL]
addendum_georg['locTokens'] = frame.locTokens

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
