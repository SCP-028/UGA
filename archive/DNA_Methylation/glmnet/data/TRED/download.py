#!python3
import requests
from bs4 import BeautifulSoup
import re

URL = "https://cb.utdallas.edu/cgi-bin/TRED/tred.cgi"
TF = {"process": "dataBrowseResult",
      "genome_id": 1,
      "browseFactor": "Browse+Factors",
      "start": 0
      }

r = requests.get(URL, params=TF)
soup = BeautifulSoup(r.text)

page_num = [num.get("value") for num in soup("option")]

TFlist = []
for num in page_num:
    TF["start"] = num
    r = requests.get(URL, params=TF)
    soup = BeautifulSoup(r.text)
    for tf in soup(href=re.compile("fid=\d*$")):
        TFlist.append(tf.get_text())


TARGET = {"process": "searchTFGene",
          "factor_organism": "any",
          "target_organism": "human",
          "sel_type": "factor_name",
          "prom_quality": 1,
          "bind_quality": 1,
          "tx_search_terms": "",
          "start": 0
          }

r = requests.get(URL, params=TARGET)
soup = BeautifulSoup(r.text)

with open("./result.csv", "w") as f:
    f.write("TF\ttarget_gene\n")
    for item in TFlist:
        gene_list = []
        TARGET["tx_search_terms"] = item
        TARGET["start"] = 0
        r = requests.get(URL, params=TARGET)
        soup = BeautifulSoup(r.text)
        for gene in soup(href=re.compile("gid=\d*$")):
                gene_list.append(gene.get_text())
        page_num = [num.get("value") for num in soup("option", value=re.compile("\d+"))]
        if len(page_num) == 1:
            pass
        else:
            for num in page_num[1:len(page_num)]:
                TARGET["start"] = num
                r = requests.get(URL, params=TARGET)
                soup = BeautifulSoup(r.text)
                for gene in soup(href=re.compile("gid=\d*$")):
                    gene_list.append(gene.get_text())
        f.write("%s\t%s\n" % (item, "; ".join(gene_list)))
