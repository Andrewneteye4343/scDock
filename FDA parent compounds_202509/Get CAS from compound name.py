#!/usr/bin/env python
# coding: utf-8

import pubchempy as pcp
import re
import csv
import os

# Variable "compounds" can be a list format or you can create it by reading your .txt file.
# Suggested format: compounds = ["compound name 1", "compound name 2", "compound name 3"]
# Suggested format: with open("Your_Path_To_TXT.txt", "r", encoding="utf-8") as f:
                        # compounds = [line.strip() for line in f]
# You need to change the definition of "compounds" for your requirement.
with open("Your_Path_To_TXT.txt", "r", encoding="utf-8") as f:
                        compounds = [line.strip() for line in f]

# Split compound names with ; symbol (due to multi-components from FDA list) and remain the unique compound
new_compounds = set()
for item in compounds:
    parts = [x.strip() for x in item.split(';')]
    new_compounds.update(parts)

new_compounds = list(new_compounds)
print("Unique compounds to process:", new_compounds)

# CSV path to save your results.
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
csv_path = os.path.join(desktop, "compound.csv")

with open(csv_path, mode="w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        "Input Name", "Input CID", "Input CAS",
        "Parent Name", "Parent CID", "Parent CAS"
    ])

    for name in new_compounds:
        print("Resolving", name)
        try:
            cids = pcp.get_cids(name, 'name')
            if not cids:
                writer.writerow([name, "'Not found", "'Not found", "'Not found", "'Not found", "'Not found"])
                continue

            cid = cids[0]
            cid_str = f"'{cid}"

            compound = pcp.Compound.from_cid(cid)
            cas = "Not found"
            for syn in compound.synonyms or []:
                if re.match(r'^\d{2,7}-\d{2}-\d$', syn):
                    cas = f"'{syn}"
                    break

            # Try to get parent compounds' PubChem CID.
            parent_cid_list = pcp.get_cids(identifier=cid, namespace='cid', domain='compound', cids_type='parent')
            if parent_cid_list:
                p_cid = parent_cid_list[0]
                p_cid_str = f"'{p_cid}"

                try:
                    parent = pcp.Compound.from_cid(p_cid)
                    p_name = parent.iupac_name or (parent.synonyms[0] if parent.synonyms else "Not found")
                    p_cas = "Not found"
                    for syn in parent.synonyms or []:
                        if re.match(r'^\d{2,7}-\d{2}-\d$', syn):
                            p_cas = f"'{syn}"
                            break
                except Exception as e:
                    p_name = "Error"
                    p_cas = "Error"
                    p_cid_str = "'Error"
            else:
                p_name = "Not found"
                p_cid_str = "'Not found"
                p_cas = "Not found"

            writer.writerow([name, cid_str, cas, p_name, p_cid_str, p_cas])

        except Exception as e:
            writer.writerow([name, "'Error", "'Error", "'Error", "'Error", "'Error"])
            print("Error with", name, ":", e)

print("Done! CSV saved at:", csv_path)