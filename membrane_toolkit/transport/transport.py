# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# Transport class


class Transport(Bulk0, Int0, Mem, IntL, BulkL):
    """
    Container class for a membrane system comprising
    bulk fluids, interfaces, and membrane

    Args:
        Bulk0 (XXX): Bulk phase on x = 0 side
        Int0 (XXX): Interface between Bulk0 and Mem
        Mem (Membrane): Membrane
        IntL (XXX): Interface between Mem and BulkL
        BulkL (XXX): Bulk phase on x=L side
    """
    pass
