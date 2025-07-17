import pandas as pd


def deduplicate_clinical(df, study_priority, exclude_cols=None, store_merged_columns=None):
    df = df.copy()
    df["priority"] = df["Study ID"].map(study_priority).fillna(999)
    df_sorted = df.sort_values(by=["Sample ID", "priority"])

    if exclude_cols is None:
        exclude_cols = []

    def merge_group(group, sample_id):
        group = group.copy()
        group["Sample ID"] = sample_id  # Inject group key manually
        group_sorted = group.sort_values("priority")
        merged = group_sorted.iloc[0].copy()
        contributing_studies = set([merged["Study ID"]])
        filled_columns = set()

        for _, row in group_sorted.iloc[1:].iterrows():
            contributing_studies.add(row["Study ID"])
            for col in group.columns:
                if col in exclude_cols or col == "priority":
                    continue
                if pd.isna(merged[col]) and not pd.isna(row[col]):
                    merged[col] = row[col]
                    filled_columns.add(col)

        merged["Study ID"] = "+".join(
            sorted(contributing_studies, key=lambda x: study_priority.get(x, 999))
        )
        if store_merged_columns:
            merged[store_merged_columns] = (
                "+".join(sorted(filled_columns)) if filled_columns else ""
            )

        return pd.DataFrame([merged])

    deduped_df = (
        df_sorted.groupby("Sample ID", group_keys=False)
        .apply(lambda g: merge_group(g, g.name), include_groups=False)
        .drop(columns=["priority"])
        .reset_index(drop=True)
    )

    cols = deduped_df.columns.tolist()
    if "Sample ID" in cols:
        cols.insert(0, cols.pop(cols.index("Sample ID")))
        deduped_df = deduped_df[cols]

    return deduped_df


def test_deduplicate_clinical():
    # Sample input
    df = pd.DataFrame(
        [
            {"Sample ID": "001", "Study ID": "study_a", "Gender": "F", "Age": 34},
            {"Sample ID": "002", "Study ID": "study_b", "Gender": None, "Age": 45},
            {"Sample ID": "002", "Study ID": "study_a", "Gender": "M", "Age": None},
            {"Sample ID": "003", "Study ID": "study_c", "Gender": "F", "Age": 50},
        ]
    )

    # Priority: study_a > study_b > study_c
    priority = {"study_a": 0, "study_b": 1, "study_c": 2}

    # Run deduplication
    deduped = deduplicate_clinical(df, priority, store_merged_columns="merged_fields")

    # Check resulting number of rows
    assert len(deduped) == 3, "Should return 3 unique Sample IDs"

    # Check Sample ID 001: no merge, study_a remains
    row_001 = deduped[deduped["Sample ID"] == "001"].iloc[0]
    assert row_001["Study ID"] == "study_a"
    assert row_001["Gender"] == "F"
    assert row_001["Age"] == 34

    # Check Sample ID 002: merge happens
    row_002 = deduped[deduped["Sample ID"] == "002"].iloc[0]
    assert row_002["Study ID"] == "study_a+study_b"
    assert row_002["Gender"] == "M"  # from study_a
    assert row_002["Age"] == 45  # from study_b

    # Check Sample ID 003: no merge
    row_003 = deduped[deduped["Sample ID"] == "003"].iloc[0]
    assert row_003["Study ID"] == "study_c"
    assert row_003["Gender"] == "F"
    assert row_003["Age"] == 50

    print("✅ All tests passed.")
    print("Deduplicated DataFrame:")
    print(deduped)


if __name__ == "__main__":
    test_deduplicate_clinical()
