#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Use Named Entity Recognition (NER) to extract location
names (tokens) from texts.
"""

import pandas as pd
import numpy as np
import spacy


NER_MODEL_DIR = snakemake.params.ner_dir
SOURCE_COLUMN = snakemake.params.source_column
LONG_STR_SEP = snakemake.params.long_str_sep
SHORT_STR_SEP = snakemake.params.short_str_sep


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
    location_tokens = series.apply(
        lambda x: _extract_location_tokens(x, model))
    return location_tokens


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
        token_strings = series.apply(lambda x: sep.join(x))
    else:
        token_strings = series.apply(lambda x: sep.join(x[-tail:]))
    token_strings.replace('', np.nan, inplace=True)
    return token_strings


# Load Darwin Core data
frame = pd.read_table(
    snakemake.input[0], usecols=['catalogNumber', SOURCE_COLUMN], dtype=str
).dropna(subset=[SOURCE_COLUMN])

# Load NER model
nlp = spacy.load(NER_MODEL_DIR)
assert 'LOC' in nlp.entity.labels

tokens = get_location_tokens(frame[SOURCE_COLUMN], nlp)
tokens_str_long = tokens_to_strings(tokens, sep=LONG_STR_SEP)
tokens_str_short = tokens_to_strings(tokens, sep=SHORT_STR_SEP, tail=3)

tokens_frame = pd.DataFrame(
    {
        'catalogNumber': frame.catalogNumber,
        'locTokens': tokens,
        'locTokensStrLong': tokens_str_long,
        'locTokensStrShort': tokens_str_short
    }, index=frame.index).dropna(subset=['locTokensStrLong'])

tokens_frame.to_csv(snakemake.output[0], index=False, sep='\t')
