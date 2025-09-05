#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os, sys, requests, time

def download_alphafold_pdb(af_id, output_dir="downloads"):
    output_dir = os.path.abspath(output_dir)  # 確保是絕對路徑
    os.makedirs(output_dir, exist_ok=True)
    url = f"https://alphafold.ebi.ac.uk/files/{af_id}-model_v4.pdb"
    output_file = os.path.join(output_dir, f"{af_id}.pdb")
    
    r = requests.get(url)
    if r.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(r.content)
        print(f"✅ Downloaded {af_id} to {output_file}")
    else:
        print(f"❌ Failed to download {af_id}, HTTP {r.status_code}")
    time.sleep(5)

if __name__ == "__main__":
    af_id, out_dir = sys.argv[1], sys.argv[2]
    download_alphafold_pdb(af_id, out_dir)

