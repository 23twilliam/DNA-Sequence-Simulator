import pytest
from src.io_functions.writers import write_fasta, write_fastq

@pytest.fixture
def sample_reads():
    return [
        {
            "id": "seq1_read_1",
            "description": "simulated_read platform=Illumina start=1 end=100",
            "sequence": "ATCGATCG",
            "qualities": [30, 30, 30, 30, 30, 30, 30, 30],
        },
        {
            "id": "seq1_read_2",
            "description": "simulated_read platform=Illumina start=50 end=150",
            "sequence": "GGCCGGCC",
            "qualities": [20, 25, 30, 35, 20, 25, 30, 35],
        }
    ]

def test_fasta_creates_file(tmp_path, sample_reads):
    out = tmp_path / "out.fasta"
    write_fasta(sample_reads, str(out))
    assert out.exists()

def test_fasta_line_count_correct(tmp_path, sample_reads):
    out = tmp_path / "out.fasta"
    write_fasta(sample_reads, str(out))
    lines = out.read_text().split("\n")
    assert len(lines) == len(sample_reads) * 2

def test_fasta_sequence_line(tmp_path, sample_reads):
    out = tmp_path / "out.fasta"
    write_fasta(sample_reads, str(out))
    lines = out.read_text().split("\n")
    assert lines[1] == "ATCGATCG"

def test_fasta_all_reads_written(tmp_path, sample_reads):
    out = tmp_path / "out.fasta"
    write_fasta(sample_reads, str(out))
    result = out.read_text()
    assert "seq1_read_1" in result
    assert "seq1_read_2" in result

########################################################################################################################

def test_fastq_creates_file(tmp_path, sample_reads):
    out = tmp_path / "out.fastq"
    write_fastq(sample_reads, str(out))
    assert out.exists()

def test_fastq_line_count_correct(tmp_path, sample_reads):
    out = tmp_path / "out.fastq"
    write_fastq(sample_reads, str(out))
    lines = out.read_text().split("\n")
    assert len(lines) == len(sample_reads) * 4

def test_fastq_header_format_correct(tmp_path, sample_reads):
    out = tmp_path / "out.fastq"
    write_fastq(sample_reads, str(out))
    lines = out.read_text().split("\n")

    assert lines[0].startswith("@")
    assert "seq1_read_1" in lines[0]
    assert "simulated_read" in lines[0]
    assert lines[4].startswith("@")
    assert "seq1_read_2" in lines[4]
    assert "simulated_read" in lines[4]

def test_fastq_sequence_line(tmp_path, sample_reads):
    out = tmp_path / "out.fastq"
    write_fastq(sample_reads, str(out))
    lines = out.read_text().split("\n")
    assert lines[1] == "ATCGATCG"

def test_fastq_plus_line(tmp_path, sample_reads):
    out = tmp_path / "out.fastq"
    write_fastq(sample_reads, str(out))
    lines = out.read_text().split("\n")
    assert lines[2] == "+"

def test_fastq_quality_encoding(tmp_path, sample_reads):
    out = tmp_path / "out.fastq"
    write_fastq(sample_reads, str(out))
    lines = out.read_text().split("\n")
    assert lines[3] == "?" * 8  # Q30 + 33 = ASCII 63 = '?' all should be Q30

def test_fastq_quality_sequence_length_matches(tmp_path, sample_reads):
    out = tmp_path / "out.fastq"
    write_fastq(sample_reads, str(out))
    lines = out.read_text().split("\n")
    for i in range(len(sample_reads)):
        seq_line = lines[i * 4 + 1]
        qual_line = lines[i * 4 + 3]
        assert len(seq_line) == len(qual_line)

def test_fastq_all_reads_written(tmp_path, sample_reads):
    out = tmp_path / "out.fastq"
    write_fastq(sample_reads, str(out))
    result = out.read_text()
    assert "seq1_read_1" in result
    assert "seq1_read_2" in result

def test_fastq_empty_reads_list(tmp_path):
    out = tmp_path / "out.fastq"
    write_fastq([], str(out))
    assert out.exists()
    assert out.read_text() == ""


def test_fastq_sequence_preserved(tmp_path, sample_reads):
    # Write reads to fastq then parse back and make sure sequence is unchanged
    out = tmp_path / "output.fastq"
    write_fastq(sample_reads, str(out))

    result = out.read_text().split("\n")
    written_sequences = []
    for i in range(len(sample_reads)):
        written_sequences.append(result[i * 4 + 1])

    assert len(written_sequences) == len(sample_reads) # if this assertion fails, the zip might still work but there has been an error
    for original, written in zip(sample_reads, written_sequences):
        assert original["sequence"] == written


def test_fastq_roundtrip_quality_preserved(tmp_path, sample_reads):
    # Write reads to fastq then convert quality back to numerical (from ASCII) and check scores are correct
    out = tmp_path / "output.fastq"
    write_fastq(sample_reads, str(out))

    lines = out.read_text().split("\n")
    for idx, read in enumerate(sample_reads):
        quality_line = lines[idx * 4 + 3]
        numerical_score = [ord(q) - 33 for q in quality_line]
        assert numerical_score == read["qualities"]