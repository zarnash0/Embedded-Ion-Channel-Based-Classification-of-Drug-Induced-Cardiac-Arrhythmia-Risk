from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings
import pandas as pd, time

Settings.Instance().MAX_LIMIT = 100
Settings.Instance().TIMEOUT = 60

act = new_client.activity
act.set_format('json')

FIELDS = ["canonical_smiles", "standard_value", "standard_units",
          "standard_relation", "molecule_chembl_id", "molecule_pref_name"]

q = act.filter(target_chembl_id="CHEMBL240", standard_type="IC50").only(FIELDS)
total = len(q)
print("Total records:", total)

rows, offset = [], 0
while offset < total:
    for attempt in range(5):
        try:
            chunk = list(q[offset:offset + 100])
            rows.extend(chunk)
            print(f"  {offset + len(chunk)} / {total}")
            break
        except Exception as e:
            print(f"  retry {attempt + 1} at offset {offset}: {type(e).__name__}")
            time.sleep(5 * (attempt + 1))
    else:
        print(f"  SKIPPING offset {offset} after 5 failures")
    offset += 100
    pd.DataFrame(rows).to_csv("herg_raw.csv", index=False)

print(f"Done. {len(rows)} records in herg_raw.csv")