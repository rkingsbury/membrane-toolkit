# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Tests for Drones
"""
import pytest
from pathlib import Path
from pandas import DataFrame

from membrane_toolkit.pipeline.drones import (
    ExptDrone,
    PermselectivityDrone,
)
from maggma.stores import MemoryStore


def test_expt_drone():
    store = MemoryStore(key="record_key")
    path = Path(__file__).absolute().parents[3] / "examples/permselectivity_data"
    config = "apparent_permselectivity.yaml"

    drone = ExptDrone(path, store, config)
    drone.run()
    assert len([e for e in drone.store.query()]) == 6


def test_permselectivity_drone():
    path = Path(__file__).absolute().parents[3] / "examples/permselectivity_data"
    drone = PermselectivityDrone(path)
    drone.run()
    assert len([e for e in drone.store.query()]) == 6
    assert isinstance(drone.store.as_df(), DataFrame)
    # drop_cols should result in "state_hash" being removed from the DataFrame
    with pytest.raises(KeyError):
        drone.store.as_df()["state_hash"]
