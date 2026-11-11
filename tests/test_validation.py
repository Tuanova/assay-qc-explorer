from assay_qc.generator import generate_dataset
from assay_qc.schema import BLANK
from assay_qc.validation import validate_data


def example():
    data, _ = generate_dataset(n_plates=2, seed=1)
    return data


def test_generated_data_is_valid_even_with_missing_signals():
    data, _ = generate_dataset(n_plates=10, seed=3, clean_fraction=0.0)
    assert data["signal"].isna().any()
    assert validate_data(data) == ([], [])


def test_missing_column_is_reported():
    errors, _ = validate_data(example().drop(columns="signal"))
    assert errors == ["Missing required columns: signal."]


def test_invalid_well_address_is_reported():
    data = example()
    data.loc[0, "well"] = "A01"
    errors, _ = validate_data(data)
    assert any("A01" in message for message in errors)


def test_duplicate_well_is_reported():
    data = example()
    data.loc[1, "well"] = data.loc[0, "well"]
    errors, _ = validate_data(data)
    assert any("more than once" in message for message in errors)


def test_text_in_signal_is_reported():
    data = example().astype({"signal": object})
    data.loc[0, "signal"] = "OVER"
    errors, _ = validate_data(data)
    assert any("non-numeric" in message for message in errors)


def test_plate_without_blank_gives_a_warning_not_an_error():
    data = example()
    data = data[data["sample_type"] != BLANK]
    errors, warnings = validate_data(data)
    assert errors == []
    assert len(warnings) == 2  # one per plate
