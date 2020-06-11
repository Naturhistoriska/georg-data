#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Use Named Entity Recognition (NER) to extract location
names (tokens) from text columns.
"""

import pandas as pd
import numpy as np
import spacy


CONFIG = snakemake.params.extraction_config


def _extract_location_tokens(text, model):
    """
    Extract location tokens from string.

    Parameters
    ----------
    text : str
    model : spaCy NER model
    """
    doc = model(text)
    location_tokens = [ent.text for ent in doc.ents if ent.label_ == 'LOC']
    return location_tokens


def get_location_tokens(series, model):
    """Return a pandas Series with extracted location tokens."""
    location_tokens = series.dropna().apply(
        lambda x: _extract_location_tokens(x, model))
    return (
        pd.Series(location_tokens, index=series.index)
        .apply(lambda d: d if isinstance(d, list) else []))


def tokens_to_strings(series, sep='', tail=None):
    """
    Return a pandas Series with concatenated tokens, optionally
    including only the last ones.

    Parameters
    ----------
    series : pandas.Series
    sep : str
        Separator to use when concatenating strings.
    tail : int
        Number of elements (from the end of the list) to include.
    """
    if tail is None:
        tokens_strings = series.apply(lambda x: sep.join(x))
    else:
        tokens_strings = series.apply(lambda x: sep.join(x[-tail:]))
    tokens_strings.replace('', np.nan, inplace=True)
    return tokens_strings


# Load Darwin Core data
source_columns = [v['source'] for _, v in CONFIG.items()]

frame = pd.read_table(
    snakemake.input[0], usecols=['occurrenceID'] + source_columns, dtype=str
).dropna(subset=source_columns, how='all')

tokens_frame = pd.DataFrame({'occurrenceID': frame.occurrenceID})

for k, v in CONFIG.items():

    nlp = spacy.load(v['nerModelDir'])  # load NER model
    assert 'LOC' in nlp.entity.labels

    tokens = get_location_tokens(frame[v['source']], nlp)
    tokens_frame[k] = tokens
    tokens_frame[k + 'StrLong'] = tokens_to_strings(
        tokens, sep=v['longStrSep'])
    tokens_frame[k + 'StrShort'] = tokens_to_strings(
        tokens, sep=v['shortStrSep'], tail=3)

tokens_frame.to_csv(snakemake.output[0], index=False, sep='\t')
