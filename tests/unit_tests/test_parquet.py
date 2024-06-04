from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import scipy

from spectrum_io.file import parquet


class TestParquet:
    """Test class to check Parquet file I/O"""

    def test_read_file(self, raw_data, tmpdir):
        output_path = Path(tmpdir / "table.parquet")
        pq.write_table(pa.Table.from_pydict(raw_data), output_path)
        df = parquet.read_file(output_path)
        pd.testing.assert_frame_equal(df, pd.DataFrame(raw_data))

    def test_write_file(self, raw_data, tmpdir):
        output_path = Path(tmpdir / "table.parquet")
        df = pd.DataFrame(raw_data)
        parquet.write_file(df, output_path)
        pd.testing.assert_frame_equal(df, pd.read_parquet(output_path))

    def test_read_write_partition(self, raw_data, tmpdir):
        """Check whether data is unmodified after being written to and then read from a partitioned dataset"""
        output_path = Path(tmpdir / "partition")
        df = pd.DataFrame(raw_data)
        parquet.write_partition([df, df], output_path, ["dataset_1", "dataset_2"])
        read_df = parquet.read_partition(output_path, "dataset_1")
        pd.testing.assert_frame_equal(read_df, df)

    def test_read_write_partition_integer_key(self, raw_data, tmpdir):
        """Check whether data is unmodified after being written to and then read from a partitioned dataset if strings of integers are used as keys, which Parquet converts to int for some reason"""
        output_path = Path(tmpdir / "partition")
        df = pd.DataFrame(raw_data)
        parquet.write_partition([df, df], output_path, ['1', '2'])
        read_df = parquet.read_partition(output_path, '1')
        pd.testing.assert_frame_equal(read_df, df)

    def test_modify_partition(self, raw_data, tmpdir):
        """Check whether file content stays the same when writing new data to the same partitioned directory"""
        output_path = Path(tmpdir / "partition")
        df = pd.DataFrame(raw_data)
        parquet.write_partition([df, df], output_path, ['1', '2'])
        parquet.write_partition([df, df, df], output_path, ['1', '2', '3'])
        read_df = parquet.read_partition(output_path, '2')
        pd.testing.assert_frame_equal(read_df, df)

    def test_write_not_implemented(self, raw_data, tmpdir):
        with pytest.raises(NotImplementedError):
            output_path = Path(tmpdir / "table.parquet")
            df = pd.DataFrame(raw_data).to_numpy()
            parquet.write_file(df, output_path)

    def test_read_write_partition_not_implemented(self, raw_data, tmpdir):
        with pytest.raises(NotImplementedError):
            output_path = Path(tmpdir / "partition")
            df = pd.DataFrame(raw_data).to_numpy()
            parquet.write_partition([df, df], output_path, ["dataset_1", "dataset_2"])

@pytest.fixture
def raw_data():
    """Simple toy MS data containing float, list[float], str, int, and list[int]"""
    return {
        "scan_number": [1, 234, 5678],
        "intensities": [
            [4e-5, 0., -1., 0., 0., -1., 0.03, 0., -1., 0.4],
            [.3, 0., -1., 1., 0., -1., 0.4, 0., -1., 0.05],
            [.04, 0., 0., 0., 0., 0., 2e-3, 0., 0., .13]
        ],
        "sequence": ["SVFLTFLR", "KTSQIFLAK", "SPVGRVTPKEWR"],
        "precursor_charge_onehot": [
            [0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
        ],
        "collision_energy_normed": [.250827308624, .288798207462, .2887064038764]
    }
