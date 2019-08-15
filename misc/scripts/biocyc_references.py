from requests_html import HTMLSession, user_agent
import pandas as pd
import re
import time
import os

df = pd.read_csv("pathway_list.csv")
urls = df.URL.tolist()
urls = [x.strip() for x in urls]
if os.path.isfile("biocyc_pubmed.csv"):
    df = pd.read_csv("biocyc_pubmed.csv")
    finished_urls = df["Pathway"].tolist()
    finished_urls = ["https://biocyc.org/HUMAN/NEW-IMAGE?type=PATHWAY&object=" + x for x in finished_urls]
    urls = [x for x in urls if x not in finished_urls]

ss = HTMLSession()
res = []
pathways = []
with open("biocyc_pubmed.csv", "a") as f:
    for pw in urls:
        headers = {
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Encoding": "gzip, deflate, br",
            "Host": "biocyc.org",
            "Accept-Language": "en-US",
            "Cache-Control": "no-cache",
            "User-Agent": user_agent()
        }
        r = ss.get(pw, headers=headers, timeout=10)
        if r.ok:
            r.html.render(wait=1, retries=3)
            # Before the h3 element that contains "Other References", get all <p> elements with class "ecoparagraph reference"
            # Get child links and their href attributes
            refs = r.html.xpath("//h3[contains(text(), 'Other References')]/preceding-sibling::p[contains(@class, 'ecoparagraph reference')]//a/@href")
            refs = [re.search(r".*pubmed/(\d+)$", x) for x in refs]
            refs = [x.group(1) for x in refs if x is not None]
            refs = ",".join(refs)
            if len(refs) > 0:
                print(f"Link: {pw}\nPubmed ID: {refs}")
                pathways.append(re.sub(r"^.*object=([^\s]*)$", r"\1", pw))
                res.append(refs)
                line = re.sub(r'^.*object=([^\s]*)$', r'\1', pw) + f',"{refs}"\n'
                f.write(line)
            else:
                f.write(f"{pw},")
        time.sleep(1)
