import pandas as pd
import pytest
from mbf_genomics.util import read_pandas
from pandas.testing import assert_frame_equal


def test_read_pandas_csv_in_xls(new_pipegraph):
    df = pd.DataFrame({"x": ["a", "b", "c"], "y": [1, 2.5, 3]})
    df.to_excel("shu.xls", index=False)
    assert_frame_equal(df, read_pandas("shu.xls"))
    df.to_csv("shu.xls", sep="\t", index=False)
    assert_frame_equal(df, read_pandas("shu.xls"))
    df.to_csv("shu.tsv", sep="\t", index=False)
    assert_frame_equal(df, read_pandas("shu.tsv"))
    df.to_csv("shu.csv", index=False)
    assert_frame_equal(df, read_pandas("shu.csv"))
    df.to_csv("shu.something", index=False)
    with pytest.raises(ValueError):
        read_pandas("shu.something")
