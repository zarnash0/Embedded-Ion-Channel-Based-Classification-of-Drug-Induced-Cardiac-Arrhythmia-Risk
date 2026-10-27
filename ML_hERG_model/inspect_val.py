import pandas as pd

SRC = r"C:\Users\User\Downloads\CiPA-Model-Validation-2018\CiPA-Model-Validation-2018\Hill_Fitting\data\manual_validationdrug_block.csv"

df = pd.read_csv(SRC)
print("Shape:", df.shape)
print("Columns:", df.columns.tolist())
print()
print(df.head(12).to_string())
print()
print("Drugs:", sorted(df.iloc[:, 0].unique()))
print()
for c in df.columns:
    if df[c].dtype == object:
        print(f"{c}: {sorted(df[c].unique())[:10]}")
    else:
        print(f"{c}: min={df[c].min()}  max={df[c].max()}")