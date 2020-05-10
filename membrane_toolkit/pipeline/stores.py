from maggma.stores import MemoryStore
from pandas import DataFrame


class PandasStore(MemoryStore):
    """
    Store that creates pandas dataframes from parsed data. At present this is just a
    lightweighte extension of MemoryStore.
    """
    def __init__(self, key: str):
        """
        Args:
            key: record key.
        """
        super().__init__(key=key)

    def as_df(self, drop_cols=["_id", "state_hash"]):
        """
        Return a Pandas dataframe representation of the data parsed by the Drone.

        Args:
            drop_cols: [str] List of column names to drop from the dataframe before returning.
                Default: ["_id", "state_hash"]
        """
        df = DataFrame(list(self.query()))
        df = df.drop(drop_cols, axis=1)
        return df
