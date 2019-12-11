"""
Placeholder for methods that facilitate calculating membrane properties from
experimental data. e.g. processing volume change vs. time into permeability.
"""
from abc import ABC
from monty.json import MSONable


class ExptBase(ABC, MSONable):
    """
    Abstract base class for analyzing experimental data. All
    derived classes must implement ...

    In general, these should read data from a .csv or .xlsx file
    (that's one method)

    Experimental parameters (pressure, geometry, etc.) need
    to be input somehow. That's probably another method or
    property of the class

    Calculation methods / options (e.g., whether to activity
    correct or use frame-of-reference corrections) have to be specified

    the process_data actually carries out the calculation and
    creates some type of output structure or object.

    The whole thing should be well-documented (perhaps via an explain dict)
    and serializable.

    """

    @abstractmethod
    def process_data(self):
        """
        Process data into membrane characteristics.
        """
        pass

    @abstractmethod
    def load_data(self):
        """
        Read experimental data from a file.
        """

    @abstractmethod
    def setup_test(self):
        """
        
        """
        pass

    @classmethod
    def from_dict(cls, d) -> 'ExptBase':
        """
        :param d: Dict representation.
        :return: ExptBase
        """
        # dec = MontyDecoder()
        # return cls(d["composition"], d["energy"], d["correction"],
        #            parameters={k: dec.process_decoded(v)
        #                        for k, v in d.get("parameters", {}).items()},
        #            data={k: dec.process_decoded(v)
        #                  for k, v in d.get("data", {}).items()},
        #            entry_id=d.get("entry_id", None))
        pass

    def as_dict(self) -> dict:
        """
        :return: MSONable dict.
        """
        pass


class DeadEndExpt(ExptBase):
    """
    Class for processing dead-end cell data.
    """


class CrossflowExpt(ExptBase):
    """
    Class for processing crossflow data.
    """


class DiffusionCellExpt(ExptBase):
    """
    Class for processing diffusion cell data.
    """
