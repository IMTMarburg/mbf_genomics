from mbf_genomics.annotator import Annotator, FromFile
from typing import Dict, List
from pypipegraph import Job
from mbf_genomics import DelayedDataFrame
from pandas import DataFrame
import pandas as pd
import hashlib


class Description(Annotator):
    """Add the description for the genes from genome.

    @genome may be None (default), then the ddf is queried for a '.genome'
    Requires a genome with df_genes_meta - e.g. EnsemblGenomes
    """

    columns = ["description"]

    def __init__(self, genome=None):
        self.genome = genome

    def calc_ddf(self, ddf):
        if self.genome is None:
            try:
                genome = ddf.genome
            except AttributeError:
                raise AttributeError(
                    "ddf had no .genome and no genome was passed to Description"
                )
        else:
            genome = self.genome
        lookup = dict(genome.df_genes_meta["description"].items())
        result = []
        for gene_stable_id in ddf.df["gene_stable_id"]:
            result.append(lookup.get(gene_stable_id, ""))
        return pd.Series(result, index=ddf.df.index)


def GeneStrandedSalmon(*args, **kwargs):
    """Deprecated. use anno_tag_counts.Salmon
    """
    raise NotImplementedError("Deprecated. Use anno_tag_counts.Salmon")


# FromFile forwarded to mbf_genomics.annotator.FromFile
FromFile = FromFile


class TMM(Annotator):
    """
    Calculates the TMM normalization from edgeR on some raw counts.
    
    Parameters
    ----------
    raw : Dict[str, Annotator]
        Dictionary of raw count annotator for all samples.
    dependencies : List[Job], optional
        List of additional dependencies, by default [].
    samples_to_group : Dict[str, str], optional
        A dictionary sample name to group name, by default None.
    """

    def __init__(
        self,
        raw: Dict[str, Annotator],
        dependencies: List[Job] = None,
        samples_to_group: Dict[str, str] = None,
    ):
        """Constructor."""
        import mbf_r # so we fail early

        self.sample_column_lookup = {}
        for sample_name in raw:
            self.sample_column_lookup[
                sample_name
            ] = f"{raw[sample_name].columns[0]} TMM"
        self.columns = list(self.sample_column_lookup.values())
        self.dependencies = []
        if dependencies is not None:
            self.dependencies = dependencies
        self.raw = raw
        self.samples_to_group = samples_to_group
        self.cache_name = hashlib.md5(self.columns[0].encode("utf-8")).hexdigest()

    def calc_ddf(self, ddf: DelayedDataFrame) -> DataFrame:
        """
        Calculates TMM columns to be added to the ddf instance.

        TMM columns are calculated using edgeR with all samples given in self.raw.

        Parameters
        ----------
        ddf : DelayedDataFrame
            The DelayedDataFrame instance to be annotated.

        Returns
        -------
        DataFrame
            A dataframe containing TMM normalized columns for each
        """
        raw_columns = [self.raw[sample_name].columns[0] for sample_name in self.raw]
        df = ddf.df[raw_columns]
        df_res = self.call_edgeR(df)
        rename = {}
        for col in df_res.columns:
            rename[col] = f"{col} TMM"
        df_res = df_res.rename(columns=rename)
        return df_res

    def call_edgeR(self, df_counts: DataFrame) -> DataFrame:
        """
        Call to edgeR via r2py to get TMM (trimmed mean of M-values)
        normalization for raw counts.

        Prepare the edgeR input in python and call edgeR calcNormFactors via
        r2py. The TMM normalized values are returned in a DataFrame which
        is converted back to pandas DataFrame via r2py.

        Parameters
        ----------
        df_counts : DataFrame
            The dataframe containing the raw counts.

        Returns
        -------
        DataFrame
            A dataframe with TMM values (trimmed mean of M-values).
        """
        import mbf_r
        import rpy2.robjects as ro
        # import rpy2.robjects.numpy2ri as numpy2ri

        ro.r("library(edgeR)")
        df_input = df_counts
        to_df = {"lib.size": df_input.sum(axis=0)}
        if self.samples_to_group is not None:
            to_df["group"] = [
                self.samples_to_group[sample_name]
                for sample_name in self.samples_to_group
            ]
        df_samples = pd.DataFrame(to_df)
        r_counts = mbf_r.convert_dataframe_to_r(df_input)
        r_samples = mbf_r.convert_dataframe_to_r(df_samples)
        y = ro.r("DGEList")(counts=r_counts, samples=r_samples,)
        # apply TMM normalization
        y = ro.r("calcNormFactors")(y)
        norm_factors = ro.r("function(y){y$samples}")(y)
        result = mbf_r.convert_dataframe_from_r(norm_factors)
        df = df_counts.divide(result["norm.factors"].values)
        df = df.reset_index(drop=True)
        return df

    def deps(self, ddf) -> List[Job]:
        """Return ppg.jobs"""
        return self.dependencies

    def dep_annos(self) -> List[Annotator]:
        """Return other annotators"""
        return list(self.raw.values())
