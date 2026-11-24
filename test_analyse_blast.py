from analyse_blast import column_matches, identity_by_window, identity_by_gene


class FakeHsp:
    """A tiny alignment: 10 columns, one mismatch and one gap in the query."""
    query_start = 1
    query = "ACGTACGT-A"
    sbjct = "ACGTTCGTAA"


def test_column_matches_counts_query_positions():
    columns = column_matches(FakeHsp())
    assert len(columns) == 10
    assert columns[0] == (1, True)
    assert columns[4] == (5, False)      # A against T is a mismatch
    assert columns[8][0] == 8            # a gap does not move the query position


def test_identity_by_window():
    columns = column_matches(FakeHsp())
    rows = identity_by_window(columns, window=5)
    assert rows[0]["percent_identity"] == 80.0   # 4 of 5 columns match


def test_identity_by_gene():
    columns = column_matches(FakeHsp())
    rows = identity_by_gene(columns, genes={"test": (1, 4)})
    assert rows[0]["percent_identity"] == 100.0
