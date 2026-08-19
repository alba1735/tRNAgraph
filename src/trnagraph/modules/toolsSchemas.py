'''
Pydantic data models shared across tRNAgraph's variant-selection and size-split-merge
machinery (see toolsTG.py's parse_variant()/build_variant_view() and adataBuild.py's
compute_variant_contribution()/merge_variant_into_adata()). Kept in a dedicated module,
following the tools*.py convention, so these schemas have one canonical, validated definition
instead of being duplicated as plain dataclasses in the modules that use them.
'''

import re
from typing import Any, Dict, Optional, Literal
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_validator

_TAG_PATTERN = re.compile(r'^(full|[uo]\d+)$')


class VariantTag(BaseModel):
    '''
    Parsed representation of a `--variant` CLI string, e.g. "raw:u60" -> norm='raw', tag='u60'.
    '''
    model_config = ConfigDict(frozen=True)

    raw: str
    norm: Literal['norm', 'raw', 'allfeatures', 'vst']
    tag: str

    @field_validator('tag')
    @classmethod
    def _validate_tag_format(cls, v: str) -> str:
        if not _TAG_PATTERN.match(v):
            raise ValueError(f"Split tag '{v}' must be 'full' or match '<u|o><cutoff>' (e.g. 'u60', 'o50').")
        return v


class VariantContribution(BaseModel):
    '''
    Coverage matrices, numeric per-obs columns, and count/sizefactor summaries for one
    size-split variant, computed by AnnDataBuilder.compute_variant_contribution() and merged
    into an existing target AnnData object by adataBuild.merge_variant_into_adata().
    '''
    model_config = ConfigDict(arbitrary_types_allowed=True)

    x_raw: pd.DataFrame
    x_norm: pd.DataFrame
    x_norm_allfeatures: pd.DataFrame
    x_vst: Optional[np.ndarray] = None
    obsm_counts: pd.DataFrame
    sizefactors_trna: Dict[str, Any]
    sizefactors_allfeatures: Dict[str, Any]
    type_counts: pd.DataFrame
    type_real_counts: pd.DataFrame
    amino_counts: pd.DataFrame
    anticodon_counts: pd.DataFrame
    nontrna_counts: pd.DataFrame
