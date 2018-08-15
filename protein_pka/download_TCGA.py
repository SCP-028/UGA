import glob
import pandas as pd


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


df = pd.read_table("/data/TCGA/manifest/gdc_manifest_count.2018-08-14.txt")
url = "https://api.gdc.cancer.gov/data/"
urls = [url + ",".join(uuid) for uuid in chunks(df["id"].unique(), 100)]
with open("/data/TCGA/count.urls", "w") as f:
    f.write("\n".join(urls))


df = pd.read_table("/data/TCGA/manifest/gdc_manifest_FPKM.2018-08-14.txt")
url = "https://api.gdc.cancer.gov/data/"
urls = [url + ",".join(uuid) for uuid in chunks(df["id"].unique(), 100)]
with open("/data/TCGA/FPKM.urls", "w") as f:
    f.write("\n".join(urls))

"""bash
mkdir count && cd count
mv ../count.urls ./
aria2c -i count.urls
for f in ./*.gz; do tar xzf $f; done
for f in ./*/*.gz; do gunzip $f; done

cd ..
mkdir FPKM && cd FPKM
mv ../FPKM.urls ./
aria2c -i FPKM.urls
for f in ./*.gz; do tar xzf $f; done
for f in ./*/*.gz; do gunzip $f; done
"""

filenames = glob.glob("/data/TCGA/count/*/*.counts")
f = filenames[0]
df = pd.read_table(f, header=None, index_col=0)
ans = pd.DataFrame(index=df.index)
for f in filenames:
    df = pd.read_table(f, header=None, index_col=0)
    df.columns = [f[17:-13]]
    ans = pd.concat([ans, df], axis=1)
