from abc import ABC
import pandas as pd
import pypipegraph as ppg

annotator_singletons = {}


class Annotator(ABC):
    def __new__(cls, *args, **kwargs):
        cn = cls.__name__
        if ppg.util.global_pipegraph:
            if not hasattr(ppg.util.global_pipegraph, "_annotator_singleton_dict"):
                ppg.util.global_pipegraph._annotator_singleton_dict = {}
            singleton_dict = ppg.util.global_pipegraph._annotator_singleton_dict
        else:
            singleton_dict = annotator_singletons
        if not cn in singleton_dict:
            singleton_dict[cn] = {}
        key = {"args": args}
        key.update(kwargs)
        key = tuple(sorted(key.items()))
        if not key in singleton_dict[cn]:
            singleton_dict[cn][key] = object.__new__(cls)
        return singleton_dict[cn][key]

    def __hash__(self):
        return hash(self.get_cache_name())

    def get_cache_name(self):
        if not isinstance(self.columns, list):
            raise ValueError("Columns was not a list")
        if hasattr(self, "cache_name"):
            return self.cache_name
        else:
            return self.columns[0]

    def calc(self, df):
        raise NotImplementedError()

    def deps(self):
        """Return ppg.jobs"""
        return []

    def dep_annos(self):
        """Return other annotators"""
        return []


class Constant(Annotator):
    def __init__(self, column_name, value):
        self.columns = [column_name]
        self.value = value

    def calc(self, df):
        return pd.DataFrame({self.columns[0]: self.value}, index=df.index)
