import pytest
from src.io_functions.readers import read_sequences

def test_single_sequence(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nATCGATCG\n")
    result = read_sequences(str(fasta))
    assert result[">seq1"] == "ATCGATCG"

def test_multiple_sequences(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nATCGATCG\n>seq2\nATCGATCC\n")
    result = read_sequences(str(fasta))
    assert len(result) == 2
    assert result[">seq1"] == "ATCGATCG"
    assert result[">seq2"] == "ATCGATCC"

def test_multiline_sequence(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nATCGATCG\nATCGATCC\n")
    result = read_sequences(str(fasta))
    assert result[">seq1"] == "ATCGATCGATCGATCC"

def test_header_saved_properly(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1 [organism=Homo sapiens] [GeneID=1] [chromosome=1]\nATCG\n")
    result = read_sequences(str(fasta))
    assert ">seq1 [organism=Homo sapiens] [GeneID=1] [chromosome=1]" in result

def test_read_returns_dict(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nATCGATCG\n")
    result = read_sequences(str(fasta))
    assert isinstance(result, dict)

def test_whitespace_stripped_from_sequence(tmp_path):
    # Lines with whitespace should be handled
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nATCG   \n  ATCG\n")
    result = read_sequences(str(fasta))
    assert "ATCG" in result[">seq1"]

def test_sequence_with_n_bases(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">seq1\nATCGNNNNATCG\n")
    result = read_sequences(str(fasta))
    assert result[">seq1"] == "ATCGNNNNATCG"

def test_long_sequence(tmp_path):
    fasta = tmp_path / "test.fasta"
    seq = "ATCG" * 10000
    fasta.write_text(f">seq1\n{seq}\n")
    result = read_sequences(str(fasta))
    assert len(result[">seq1"]) == 40000

def test_file_not_found_raises():
    with pytest.raises((FileNotFoundError, OSError)):
        read_sequences("nonexistent_file.fasta")