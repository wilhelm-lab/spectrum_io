import shutil
import sys
import tempfile
import unittest
from contextlib import nullcontext
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy

if "typeguard" in sys.modules:
    from typeguard import suppress_type_checks

from spectrum_io.file import parquet


class TestParquet(unittest.TestCase):
    """Test class to check Parquet file I/O."""

    def setUp(self):  # noqa: D102
        # Simple toy MS data containing float, list[float], str, int, and list[int]
        self.raw_data = {
            "scan_number": [1, 234, 5678],
            "intensities": [
                [4e-5, 0.0, -1.0, 0.0, 0.0, -1.0, 0.03, 0.0, -1.0, 0.4],
                [0.3, 0.0, -1.0, 1.0, 0.0, -1.0, 0.4, 0.0, -1.0, 0.05],
                [0.04, 0.0, 0.0, 0.0, 0.0, 0.0, 2e-3, 0.0, 0.0, 0.13],
            ],
            "sequence": ["SVFLTFLR", "KTSQIFLAK", "SPVGRVTPKEWR"],
            "precursor_charge_onehot": [
                [0, 1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
            ],
            "collision_energy_normed": [0.250827308624, 0.288798207462, 0.2887064038764],
        }
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):  # noqa: D102
        shutil.rmtree(self.temp_dir)

    def test_read_file(self):
        """Test read operation for a single dataset."""
        output_path = self.temp_dir / "table.parquet"
        pq.write_table(pa.Table.from_pydict(self.raw_data), output_path)
        df = parquet.read_file(output_path)
        pd.testing.assert_frame_equal(df, pd.DataFrame(self.raw_data))

    def test_write_file(self):
        """Check write operation for a single dataset."""
        output_path = self.temp_dir / "table.parquet"
        df = pd.DataFrame(self.raw_data)
        parquet.write_file(df, output_path)
        pd.testing.assert_frame_equal(df, pd.read_parquet(output_path))

    def test_read_write_partition(self):
        """Check whether data is unmodified after being written to and then read from a partitioned dataset."""
        output_path = self.temp_dir / "partition"
        df = pd.DataFrame(self.raw_data)
        parquet.write_partition([df, df], output_path, ["dataset_1", "dataset_2"])
        read_df = parquet.read_partition(output_path, "dataset_1")
        pd.testing.assert_frame_equal(read_df, df)

    def test_read_write_partition_integer_key(self):
        """Check whether Parquet's under-the-hood conversion of string to integer keys is handled seamlessly."""
        output_path = self.temp_dir / "partition"
        df = pd.DataFrame(self.raw_data)
        parquet.write_partition([df, df], output_path, ["1", "2"])
        read_df = parquet.read_partition(output_path, "1")
        pd.testing.assert_frame_equal(read_df, df)

    def test_modify_partition(self):
        """Check whether file content stays the same when writing new data to the same partitioned directory."""
        output_path = self.temp_dir / "partition"
        df = pd.DataFrame(self.raw_data)
        parquet.write_partition([df, df], output_path, ["1", "2"])
        parquet.write_partition([df, df, df], output_path, ["1", "2", "3"])
        read_df = parquet.read_partition(output_path, "2")
        pd.testing.assert_frame_equal(read_df, df)
