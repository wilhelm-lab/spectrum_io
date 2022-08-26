import pytest
import numpy as np
import pandas as pd

import prosit_io.search_result.maxquant as mq


class TestAddTMTMod:
    def test_add_tmt_mod(self):
        assert mq.MaxQuant.add_tmt_mod(1.0, '[UNIMOD:2016]ABC[UNIMOD:4]K[UNIMOD:2016]', '[UNIMOD:2016]') == 1.0 + 2*304.207146
